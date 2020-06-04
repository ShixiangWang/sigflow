#!/usr/bin/env Rscript

'=================================================================
sigflow: Streamline Analysis Workflow for Mutational Signatures.

Author: Shixiang Wang (wangshx@shanghaitech.edu.cn)
Copyright: MIT@2020

Usage:
  sigflow.R extract --input=<input> --output=<output>
  sigflow.R (-h | --help)
  sigflow.R --version

Options:
  -h --help     Show help message.
  --version     Show version.
  -i <input>, --input <file>  input file path.
  -o <output>, --output <file>  output directory path.

=================================================================
' -> doc

library(docopt)
arguments <- docopt(doc, version = 'sigflow v0.1\n')
print(arguments)
