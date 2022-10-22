def run_script(script, stdin=None):
  """Returns (stdout, stderr), raises error on non-zero return code"""
  import subprocess
  # Note: by using a list here (['bash', ...]) you avoid quoting issues, as the 
  # arguments are passed in exactly this order (spaces, quotes, and newlines won't
  # cause problems):
  proc = subprocess.Popen(['bash', '-c', script],
    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    stdin=subprocess.PIPE)
  stdout, stderr = proc.communicate()
  if proc.returncode:
    raise ScriptException(proc.returncode, stdout, stderr, script)
  return stdout, stderr

class ScriptException(Exception):
  def __init__(self, returncode, stdout, stderr, script):
    self.returncode = returncode
    self.stdout = stdout
    self.stderr = stderr
    Exception().__init__('Error in script')

name_pref = "ekrt_like"
N_events = 5000
K_factors = [4]
M_factors = [0.1, 0.2, 0.3, 0.4, 0.5, 0.65, 0.8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20]
factors = [(k,m) for k in K_factors for m in M_factors]
with open("params_template", "r") as template_f:
  template = template_f.read()
  for (K, M) in factors:
    full_name = name_pref+"_K="+str(K)+"_M="+str(M)
    copy_t = template
    with open(full_name, "w") as setting_f:
      print(copy_t.format(
      full_name, 
      "false", 
      "false", 
      "true", 
      "false", 
      "false", 
      "false", 
      "false", 
      "false", 
      "true", 
      "false", 
      "true", 
      str(N_events),
      str(20),
      str(1.0),
      str(K),
      str(M)),file=setting_f)
    run_script("./bin/mcaa "+full_name)
    full_name += "_nPDF"
    copy_t = template
    with open(full_name, "w") as setting_f:
      print(copy_t.format(
      full_name, 
      "false", 
      "false", 
      "true", 
      "true", 
      "false", 
      "false", 
      "false", 
      "false", 
      "true", 
      "false", 
      "true", 
      str(N_events),
      str(20),
      str(1.0),
      str(K),
      str(M)),file=setting_f)
    run_script("./bin/mcaa "+full_name)