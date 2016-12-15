# README-wrappers

This document describes design principles for writing wrappers to third-party tools, for internal code consistency and flexibility.

## Introduction

The `pycits` pipeline relies almost entirely on third-party tools. For future flexibility and ease of coding, it is useful to have a single core style for wrappers, so that the cognitive load of using a range of tools is minimised. This makes development much easier, and conveniently abstracts out common features/functions of wrappers (and command-line tools in general) into reusable code.

## General design

### Wrapper classes

* Wrapper classes should be named for the tool (see below). For example, the `Swarm` tool should have a module named `swarm.py` that contains a `Swarm()` class, through which programmatic interactions with the `Swarm` tool will be mediated.
* The wrapper class should provide a `.run()` method, which can be used to execute a single instance of the command-line tool.
* The `.run()` method will take a `dry_run` argument (default: `False`) which, when set to `True`, causes the method to return the command-line that would have been run, but not to run the tool itself. This is useful for compiling iterables of commands for later submission to a scheduler such as [SGE](https://en.wikipedia.org/wiki/Oracle_Grid_Engine) or [Univa Grid Engine](https://en.wikipedia.org/wiki/Univa). **NOTE: this design decision may change - it could be clearer to make the `dry_run` a separate method**
* The `.run()` method should take a required `parameters` argument that expects a `namedtuple` `Parameters()` object (see below), defined in the same module. This encapsulates the user-modifiable parameters for the tool.
* Additional tool parameters that are currently considered conceptually distinct from the `parameters` argument are input and output file paths. **NOTE: this design decision may change - separation of file paths and parameters is flexible, but it could be confusing and arbitrary**
* The `.run()` method will return a `Results()` object (see below), which is defined in the same module.
* The class will provide a `.__build_cmd()` private method that constructs the command-line to be run, and places it in the `self._cmd` class attribute, where it can be returned (in a dry run), executed (e.g. with `subprocess.run()`) or returned as part of the `Results()` object (see below).

### Parameters class

* Each module should contain a `Parameters()` class, which is a `namedtuple`. This defines a set of user-modifiable parameters, with default values. 

By having a consistent `Parameters()` object, we can reuse/simplify code to target a generic object, with minor modifications for specific tools, when required. This reduces cognitive load in maintenance and design of wrappers and integration code for new tools. The `namedtuple` pattern allows us to identify parameters directly with meaningful names, improving readbility.

### Results class

* Each module should contain a `Results()` class, which is a `namedtuple`. This defines a set of outputs expected from the tool.
* Three outputs should always be available in the `Results()` class:
  * `command`: the command-line that ran the tool (results should be repeatable by copy-pasting this into the terminal
  * `stdout`: the captured `STDOUT` from running the tool
  * `stderr`: the captured `STDERR` from running the tool
* `STDERR` and `STDOUT` are captured for two main reasons: debugging/diagnostics; and capturing useful output from some packages that only write to `STDOUT`, rather than an output file. 
* The remaining outputs may vary by tool, but are likely to include a path to an output file. **NOTE: it may be useful to standardise the names of these outputs across tools**

By having a consistent `Results()` object with predictable items returned from each tool, we can reuse/simplify code to target a generic object, with minor modifications for specific tools, when required. This reduces cognitive load in maintenance and design of wrappers and integration code for new tools. The `namedtuple` pattern allows us to identify returned values and data directly with meaningful names, improving readbility.

### Exceptions

* Each module should define an Exception that is specific for the module, and named for the module. For example, the `swarm.py` module would define the `SwarmException()` error.
* The module exception should be raised when errors are caught in this module.

Using named errors in this way makes debugging easier, as we can identify the modules from which errors derive directly. The ability to [chain exceptions](https://blog.ionelmc.ro/2014/08/03/the-most-underrated-feature-in-python-3/) in Python3 means that we lose no detail by using this model, and gain only extra information about where the error was caught.


### Location of wrappers

* Wrapper code should be contained in the `pycits` module.
* Each tool should have a module that is named for that tool, in [snake-case](https://en.wikipedia.org/wiki/Snake_case). For example, the `Trimmomatic` wrapper should be contained in the module with filename `trimmomatic.py`
* The wrapper should comprise a single object, named for the tool, in that module. For example, the `Swarm` tool should have a module named `swarm.py` that contains a `Swarm()` class, through which programmatic interactions with the `Swarm` tool will be mediated.

By having a consistent location and format for wrappers and modules like this, we reduce the cognitive load of code maintenance, and ease of modification by other users.

## Example code

For some third party tool `FooBar`, we would create - at a minimum - a module `foobar.py`, containing the `Results()`, `Parameters()`, `FoobarException()`, and `Foobar()` classes, as follows:

```python
#!/usr/bin/env python3.5
# -*- encoding: utf-8 -*-

"""
Docstring
"""

# IMPORTS HERE

# Factory class (named tuple) for FooBar returned values
Results = namedtuple("Results",
                     "command stdout stderr outfilename excludeddata")
                     
# Factory class (named tuple) for FooBar user parameters
Parameters = namedtuple("Parameters", "verbose threads iterations")
Parameters.__new__.__defaults__ = (False, 1, 100)


# Exception for this module
class FoobarError(Exception):
    """Exception raised when Foobar fails"""
    def __init__(self, msg):
    	self.message = msg
    	

# Foobar class
class Foobar(object):
    """Wrapper for FooBar"""
    def __init__(self, exe_path):
        self._exe_path = exe_path
        
    def run(self, infname, outdir, parameters, dry_run=False):
        """Run FooBar"""
        self.__build_cmd(infname, outdir, parameters)
        # builds self._cmd, and other attributes to be returned
        if dry_run:
            return(self._cmd)
        pipe = subprocess.run(self._cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=True)
        results = Results(self._cmd, pipe.stdout, pipe.stderr,
                          self._outfname, self._excluded)
        return results
        
    def __build_cmd(self, infname, outdir, parameters):
        """Construct command-line for FooBar"""
        self._outfname = os.path.join(outdir, "foobar.out")
        self._excluded = os.path.join(outdir, "excluded.tab")
        params = ["--{0} {1}".format(k, v) for k, v in
                  parameters._asdict().items() if v is not None]
        cmd = ["foobar", *params, "-o {0}".format(self._outfname),
               "-e {0}".format(self._excluded)]
        self._cmd = ' '.join(cmd)
```

## Usage

By keeping close to the design principles above, the usage of all third-party tools becomes very consistent, which should again help reduce cognitive load for maintenance, and integration of new tools.

As an example of usage, we continue with the `FooBar` tool example, illustrating the three steps of use, with the `foobar` module:

* Instantiate the tool wrapper object (`foobar.Foobar()`)
* Create a parameter set (using `foobar.Parameters()`)
* Run the tool and collect results (`.run()`)

```python
# Import module
from pycits import foobar

# Create output directory if it doesn't exist
os.makedirs("outdir", exist_ok=True)

# 1) Instantiate tool wrapper, with location of executable
fb_exe = foobar.Foobar("/usr/bin/foobar")

# 2) Create parameter set
fb_parameters = foobar.Parameters(verbose=False, iterations=50)

# 3) Run tool and collect results
fb_results = fb_exe.run("foobar_input.fas", "outdir", fb_parameters)
```

The `Results()` object (here, `fb_results`) can then be queried by keyword, e.g.

```python
>>> fb_results['command']
... "foobar --verbose False --threads 4 --iterations 50 " \
    "-o outdir/foobar.out -e outdir/excluded.tab"
```

If you wanted to create a number of alternative parameter sets to generate command-lines to be run on SGE, you could use code like the following:

```python
fb_cmds = []
for iterations in range(1000):
	params = foobar.Parameters(iterations=iterations)
	fb_exe = foobar.Foobar("foobar")
	fb_cmds.append(fb_exe.run("foobar_input.fas",
	                          "outdir_{0}".format(iterations),
	                          params, dry_run=True)
```

which would generate 1000 command lines in `fb_cmds` that could be organised to be submitted to a job scheduler.
