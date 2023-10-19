import io, os
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
from tqdm import tqdm
import random
import argparse
import subprocess
import numpy as np
import tensorflow as tf
import multiprocessing as mp
import tensorflow_probability as tfp
import warnings, logging

tf.get_logger().setLevel(logging.ERROR)

warnings.filterwarnings('ignore')
tfb = tfp.bijectors
tfd = tfp.distributions

try:
  print("Trying to find GPU...")
  if tf.test.gpu_device_name() != '/device:GPU:0':
    print('WARNING: GPU device not found.')
  else:
    print('SUCCESS: Found GPU: {}'.format(tf.test.gpu_device_name()))
except:
  pass

def joint_log_prob_4su_vec(pi_g, conversions, t_content, p_c, p_e):

    nreads = len(t_content)

    mix_prob = [pi_g, (1-pi_g)]
    mix_probs = tf.tile([mix_prob], multiples = (nreads, 1))
    rv_assignments = tfd.Categorical(probs = mix_probs)

    pcpe = tf.stack([p_c, p_e])
    repprobs = tf.tile([[p_c, p_e]], multiples = (nreads, 1))
    t_content = t_content.reshape(nreads, 1)
    t_content = tf.cast(t_content, dtype = 'float32')
    binoms = tfd.Binomial(total_count = t_content, probs = repprobs)

    rv_observations = tfd.MixtureSameFamily(
        mixture_distribution = rv_assignments,
        components_distribution = binoms)

    return(
        tf.reduce_sum(rv_observations.log_prob(conversions))
    )

def subsample(incont, inconv, i=5000):
  cont = []
  conv = []
  p = i / float(len(incont))
  for a,b in zip(incont, inconv):
    if random.random() < p:
      cont.append(a)
      conv.append(b)
  return np.array(cont), conv

@tf.function(jit_compile = True)
def sample_kern(num_steps = 5000, num_burnin_steps = 1000, current_state = None, kernel = None):
  return tfp.mcmc.sample_chain(
    num_results = num_steps,
    num_burnin_steps = num_burnin_steps,
    current_state = current_state,
    kernel = kernel,
    trace_fn = None)

def hmc_sampler(log_prob_fn, bijector, line, init_num, outfolder):
    
    # Data is passed as line-by-line from the combined infile
    data = line.rstrip().decode().split('\t')
    bc, gene = data[2:4]
    p_c, p_e, numreads = data[6:9]
    p_c = float(p_c)
    p_e = float(p_e)
    numreads = int(numreads)
    
    # Changing p_c to conversion probability plus background probability
    p_c = p_c + p_e
    name = f"{bc}__{gene}"
    conversions = [int(i) for i in data[4].split(',')]
    numconvs    = np.array(conversions).sum()
    t_content   = np.array(data[5].split(','), dtype = int)

    # If no conversions are found, the pi_g estimate is set to zero
    if numconvs == 0:
        with open(os.path.join(outfolder, 'pi_g_results/%s_pig.txt' % name), 'w') as outfh:
            outfh.write('%s\t%.5f\t%.5f\t%.5f\t%d\n' % (name, 0, 0, 0, numreads))
            return
    
    # Subsample highly covered genes to 5000 reads to speed up processing
    if len(conversions) > 5000:
       t_content, conversions = subsample(t_content, conversions, i=5000)

    init_state = tf.constant(init_num, dtype = 'float32')
    unp = lambda pi_g: log_prob_fn(pi_g, conversions, t_content, p_c, p_e)

    # Set up MCMC kernel
    hmc_kernel = tfp.mcmc.HamiltonianMonteCarlo(
        target_log_prob_fn = unp,
        num_leapfrog_steps = 10,
        step_size = 0.05)

    transformed_kernel = tfp.mcmc.TransformedTransitionKernel(
        inner_kernel = hmc_kernel,
        bijector = bijector)

    kernel = tfp.mcmc.SimpleStepSizeAdaptation(
      inner_kernel=transformed_kernel,
      num_adaptation_steps=int(0.8 * 1000), target_accept_prob=0.65)

    # Sample from the MCMC kernel
    # Number of sampling steps is hard-coded to 5000 and burn-in steps to 1000
    posterior_pig = sample_kern(num_steps = 5000, num_burnin_steps = 1000,
                                current_state = init_state, kernel = kernel)

    # Calculate the credible interval of sampling result
    lowperc  = np.percentile(posterior_pig.numpy(), 5)
    mean = posterior_pig.numpy().mean()
    highperc = np.percentile(posterior_pig.numpy(), 95)

    # Write the results to file
    with open(os.path.join(outfolder, 'pi_g_results/%s_pig.txt' % name), 'w') as outfh:
      outfh.write('%s\t%.5f\t%.5f\t%.5f\t%d\n' % (name, lowperc, mean, highperc, numreads))
      return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--work_dir',  help = 'Path to working directory. This folder shall include the raw data and the barcode to index mapping')
    parser.add_argument('-n', '--num_iter',  help = 'Number of iterations before restarting the process. Needed to tackle memory leakage problems', type = int, default = 5_000)
    parser.add_argument('-p', '--n_jobs',    help = 'Number of concurrent processes', type = int, default = 20)
    parser.add_argument('-r', '--readt',     help = 'Number of reads to require for inference', type = int, default = 1)
    parser.add_argument('-o', '--outfolder', help = "Folder where the pi_g_results are saved to. Default is data", type = str, default = "data")
    parser.add_argument('-s', '--start',     help = "Start index of indata for streaming data to the HMC sampler", type = int)
    parser.add_argument('-e', '--end',       help = "End index of indata for streaming data to the HMC sampler", type = int)
    parser.add_argument('--init_guess',      help = 'Initial guess for the MCMC sampler, default = 0.1', type = float, default = 0.1)
    
    args = parser.parse_args()

    # Read arguments
    start = args.start
    end   = args.end
    work_dir  = args.work_dir
    n_jobs    = args.n_jobs
    num_iter  = args.num_iter
    init_num  = args.init_guess
    outfolder = args.outfolder
    read_threshold = args.readt

    # Initiate the outfolder for the first iteration
    if not os.path.exists(os.path.join(outfolder, 'pi_g_results')):
        os.mkdir(os.path.join(outfolder, 'pi_g_results'))

    # Generate data chunk to work with from the combined index file that is tabix:ed per line
    concatenated_index_range = " ".join([str(i) for i in range(start, end)])
    print(f"Finding data to work with in the index range {start} to {end}...")
    cmd = f"tabix {os.path.join(work_dir, 'combined_results.txt.gz')} {concatenated_index_range}"
    res = subprocess.run(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    res = io.BytesIO(res.stdout)

    # Instantiate worker pool
    pool = mp.Pool(processes = n_jobs)

    # Run jobs in parallel
    print("Running MCMC sampler for data in the given range...")
    jobs = [pool.apply_async(hmc_sampler, (joint_log_prob_4su_vec, tfb.Sigmoid(), line, init_num, outfolder, )) for line in res.readlines()]
    
    results = [job.get() for job in tqdm(jobs)]

    print("Done with chunk, restarting!")
