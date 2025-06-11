using PyCall
# to load python scripts from file
pyimport("numpy")
pyimport_conda("importlib.util", "importlib")
importlib = pyimport("importlib.util")
function load_python_scripts(filename::String)
  spec = importlib.spec_from_file_location("mymodule", filename)
  mymodule = importlib.module_from_spec(spec)
  spec.loader.exec_module(mymodule)
  return mymodule
end
