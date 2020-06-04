'sigflow: Streamline Analysis Workflow for Mutational Signatures.

Author: Shixiang Wang (wangshx@shanghaitech.edu.cn)
Copyright: MIT@2020

Usage:
  sigflow.R ship new <name>...
  sigflow.R ship <name> move <x> <y> [--speed=<kn>]
  sigflow.R ship shoot <x> <y>
  sigflow.R mine (set|remove) <x> <y> [--moored | --drifting]
  sigflow.R (-h | --help)
  sigflow.R --version

Options:
  -h --help     Show this screen.
  --version     Show version.
  --speed=<kn>  Speed in knots [default: 10].
  --moored      Moored (anchored) mine.
  --drifting    Drifting mine.

' -> doc

library(docopt)
arguments <- docopt(doc, version = 'sigflow v0.1\n')
print(arguments)
