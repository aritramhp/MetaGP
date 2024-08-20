Help()
{
   me=`basename "$0"`

   # Display Help
   echo "MetaGP - Metagenomic data analysis Pipeline"
   echo "This pipeline can perform different steps in analyzing metagenomic data. It has both CPU and GPU compatible versions."
   echo
#    echo "For copyright and licensing information, please see"
#    echo "https://github.com/zjshi/gt-pro/blob/master/LICENSE"
   echo
   echo "Syntax: $me [config|pre_exec|exec_qc|exec_taxprof|exec_funcprof|-h]"
   echo "options:"
   echo "  config                   Build a configuration file with default parameters of the tools."
   echo "  pre_exec                 Execute the pre-execution step"
   echo "  exec_qc                  Execute the quality control steps"
   echo "  exec_taxprof             Execute the taxonomy profiling"
   echo "  exec_funcprof            Execute the functional profiling"
   echo "  -h [--help]              Print this help."
   echo
   echo
}

if [ "$1" = "config" ]; then
	shift
	python make_config.py $*
elif [ "$1" = "pre_exec" ]; then
	shift
	python quality_check.py $*
elif [ "$1" = "exec_qc" ]; then
	shift
	python quality_control.py $*
elif [ "$1" = "exec_taxprof" ]; then
	shift
	python taxonomic_profiling.py $*
elif [ "$1" = "exec_funcprof" ]; then
	shift
	python func_profiling.py $*
elif [ "$1" = "-h" ] || [ "$1" = "-help" ] || [ "$1" = "help" ] || [ "$1" = "--help" ]; then
	Help
else
	Help
fi