'''
Created on Dec 14, 2012

@author: mkiyer

AssemblyLine: transcriptome meta-assembly from RNA-Seq

Copyright (C) 2012 Matthew Iyer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import logging
import argparse
import os
import sys
import xml.etree.cElementTree as etree

import assemblyline.utils
_script_dir = assemblyline.utils.__path__[0]
from assemblyline.lib.librarytable import LibraryInfo

def parse_bool(s):
    firstchar = s[0].lower()
    if (firstchar == "t") or (firstchar == "y"):
        return True
    return False

def bash_log(msg, level="DEBUG"):
    return 'echo "[`date`] - %s - %s" >&2' % (level, msg)

def bash_check_retcode(msg="ERROR", level="ERROR"):
    echo_command = bash_log(msg, level)
    return 'ERRCODE=$?; if [ $ERRCODE -gt 0 ]; then %s; exit $ERRCODE; fi' % (echo_command)

def get_pbs_header(job_name,
                   num_processors=1,
                   node_processors=1,
                   node_memory=4096,
                   pbs_script_lines=None, 
                   working_dir=None, 
                   walltime=None, 
                   pmem=None, 
                   mem=None,
                   deps=None,
                   email=None,
                   stdout_filename=None,
                   stderr_filename=None):
    '''
    job_name: string name of job
    num_processors: number of processors to submit with (cannot be greater than the number of processors per node)
    node_processors: number of cores available per node
    node_memory: amount of memory per node (MB)
    pbs_script_lines: list of PBS directives to be added to the script
    working_dir: the "working directory" of the job (allows scripts to access files using relative pathnames)
    walltime: the walltime passed to qsub
    pmem: amount of memory allocated to each processor of this job in MB
    mem: amount of total memory to split across all processors (overrides pmem setting)
    deps: 'None' if no dependencies, or a python list of job ids
    email: 'None', or string containing codes 'b', 'a', or 'e' describing when to email
    stdout_filename: string filename for storing stdout
    stderr_filename: string filename for storing stderr    
    '''    
    if pbs_script_lines is None:
        pbs_script_lines = []
    if isinstance(deps, basestring):
        deps = [deps]
    # ensure number of processors can fit on node
    num_processors = min(num_processors, node_processors)
    resource_fields = ["nodes=%d:ppn=%d" % (1, num_processors)]    
    # setup memory resources
    if mem is not None:
        if mem > node_memory:
            logging.warning("Job requested more memory than node supports (%dmb > %dmb)" % (mem, node_memory))        
        mem = min(mem, node_memory)
        resource_fields.append("mem=%dmb" % (mem))
    elif pmem is not None:
        max_pmem = float(node_memory) / num_processors
        if pmem > max_pmem:
            logging.warning("Job requested more memory-per-process (pmem) than node supports (%dmb > %dmb)" % (pmem, max_pmem))
        pmem = min(pmem, max_pmem)
        resource_fields.append("pmem=%dmb" % (pmem))
    else:
        pmem = int(round(float(node_memory) / node_processors, 0))
        mem = pmem * num_processors
        resource_fields.append("mem=%dmb" % (mem))
    # set job walltime
    if walltime is not None:
        resource_fields.append("walltime=%s" % (walltime))
    # add PBS parameters
    lines = ["#PBS -N %s" % job_name,
             "#PBS -l %s" % (",".join(resource_fields))]
    if email is not None:
        lines.append("#PBS -m %s" % (email))
    if stdout_filename is None:
        stdout_filename = "/dev/null"
    lines.append("#PBS -o %s" % (stdout_filename))
    if stderr_filename is None:
        stderr_filename = "/dev/null"        
    lines.append("#PBS -e %s" % (stderr_filename))
    if deps is not None:
        lines.append("#PBS -W depend=afterok:%s" % (":".join([d for d in deps])))    
    lines.extend(pbs_script_lines)
    if working_dir is not None: 
        lines.append("cd %s" % (working_dir))
    return lines

class PluginConfig(object):
    @staticmethod
    def from_xml(xmlfile):
        tree = etree.parse(xmlfile)  
        root = tree.getroot()
        c = PluginConfig()
        # plugins
        c.plugins = {}
        for elem in root.findall("plugin"):
            name = elem.get("name")
            c.plugins[name] = elem
        # modules
        c.modules_init_script = root.findtext("modules_init_script")
        # pbs
        c.pbs = False
        c.node_mem = None
        c.node_processors = None
        c.pbs_script_lines = []
        pbs_elem = root.find("pbs")
        has_pbs = (pbs_elem.get("use") == "yes")
        if has_pbs:
            c.pbs = True
            c.node_mem = float(pbs_elem.findtext("node_mem"))
            c.node_processors = int(pbs_elem.findtext("node_processors"))
            for line_elem in pbs_elem.findall("script_line"):
                c.pbs_script_lines.append(line_elem.text)
        return c

def htseq_count_plugin(libs, config, plugin_elem, keep_tmp, dryrun):
    # setup output directories
    output_dir = plugin_elem.findtext("output_dir")
    tasks_dir = os.path.join(output_dir, "tasks")
    log_dir = os.path.join(output_dir, "logs")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(tasks_dir):
        os.makedirs(tasks_dir)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    # htseq parameters
    num_processors = int(plugin_elem.findtext("num_processors"))
    copy_bam = parse_bool(plugin_elem.findtext("copy_bam", "no"))
    run_as_pe = parse_bool(plugin_elem.findtext("pe", "no"))    
    gtf_file = plugin_elem.findtext("gtf_file")
    extra_args = []
    for arg_elem in plugin_elem.findall("arg"):
        extra_args.append(arg_elem.text)        
    # get modules to add
    modules = []
    modules_elem = plugin_elem.find("modules")
    for elem in modules_elem.findall("module"):
        modules.append(elem.text)
    # pbs defaults
    PBS_JOB_MEM = 4000
    PBS_JOB_WALLTIME = "20:00:00"
    # parse library table
    for lib in libs:
        #
        # build up sequence of commands
        #
        shell_commands = []
        shell_commands.append("#!/bin/sh")
        #
        # get pbs header
        #
        if config.pbs:
            stdout_file = os.path.join(log_dir, "%s.stdout" % (lib.library_id))
            stderr_file = os.path.join(log_dir, "%s.stderr" % (lib.library_id))
            pbs_commands = get_pbs_header(job_name=lib.library_id,
                                          num_processors=num_processors,
                                          node_processors=config.node_processors,
                                          node_memory=config.node_mem,
                                          pbs_script_lines=config.pbs_script_lines, 
                                          walltime=PBS_JOB_WALLTIME, 
                                          mem=PBS_JOB_MEM,
                                          stdout_filename=stdout_file,
                                          stderr_filename=stderr_file)
            shell_commands.extend(pbs_commands)
        #
        # intro logging messages
        #
        shell_commands.append(bash_log("AssemblyLine RNA-Seq Pipeline", "INFO"))
        shell_commands.append(bash_log("Plugin 'htseq'", "INFO"))
        shell_commands.append(bash_log("Library '%s'" % (lib.library_id), "INFO"))
        #
        # setup environment
        #
        shell_commands.append("source %s" % (config.modules_init_script))
        for name in modules:        
            basename = name.split("/")[0]
            shell_commands.append("module rm %s" % (basename))
            shell_commands.append("module add %s" % (name))
        #
        # create directories
        #
        lib_output_dir = os.path.join(output_dir, lib.library_id)
        msg = "Creating directory: %s" % (lib_output_dir)
        shell_commands.append(bash_log(msg, "INFO"))
        shell_commands.append("mkdir -p %s" % (lib_output_dir))
        shell_commands.append(bash_check_retcode(msg))  
        tmp_dir = os.path.join(lib_output_dir, "tmp")
        msg = "Creating directory: %s" % (tmp_dir)
        shell_commands.append(bash_log(msg, "INFO"))
        shell_commands.append("mkdir -p %s" % (tmp_dir))
        shell_commands.append(bash_check_retcode(msg))
        #
        # copy bam
        #
        if copy_bam:
            inp_bam_file = os.path.join(tmp_dir, os.path.basename(lib.bam_file))
            shell_commands.append("cp %s %s" % (lib.bam_file, inp_bam_file))
            shell_commands.append(bash_check_retcode("Error copying BAM file"))
        else:
            inp_bam_file = lib.bam_file   
        #
        # run htseq-count for gene expression
        #
        msg = "Counting reads across genes with htseq-count"        
        shell_commands.append(bash_log(msg, "INFO"))
        output_file = os.path.join(lib_output_dir, "htseq.txt")
        log_file = os.path.join(lib_output_dir, 'htseq.log')
        args = ["python", os.path.join(_script_dir, "run_htseq_count.py"),
                "--tmp-dir", tmp_dir]
        if run_as_pe:
            args.append("--pe")
        for arg in extra_args:
            args.append('--arg="%s"' % arg)
        args.extend([gtf_file,
                     inp_bam_file,
                     output_file,
                     '> %s 2>&1' % (log_file)])                     
        command = ' '.join(map(str, args))
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode(msg))
        #
        # cleanup
        #
        msg = "Cleaning up tmp files"
        if keep_tmp:
            shell_commands.append(bash_log("[SKIPPED] %s because '--keep-tmp' flag set" % (msg), "INFO"))
        else:
            shell_commands.append(bash_log(msg, "INFO"))
            command = "rm -rf %s" % (tmp_dir)
            shell_commands.append(command)    
            shell_commands.append(bash_check_retcode())   
        #
        # write commands 
        #
        if dryrun:
            for command in shell_commands:
                print command
        else:
            if config.pbs:
                ext = ".pbs"
            else:
                ext = ".sh"
            filename = os.path.join(tasks_dir, "%s%s" % (lib.library_id, ext))
            f = open(filename, "w")
            for command in shell_commands:
                print >>f, command
            f.close()
    return 0

PLUGINS = {'htseq': htseq_count_plugin}

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logging.info("AssemblyLine %s" % (assemblyline.__version__))
    logging.info("----------------------------------")
    # command line parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("--dryrun", dest="dryrun", action="store_true", 
                        default=False)
    parser.add_argument("--keep-tmp", dest="keep_tmp", action="store_true", 
                        default=False)
    parser.add_argument('library_table')
    parser.add_argument("config_xml_file")
    parser.add_argument("plugin_name")
    args = parser.parse_args()
    # read configuration
    logging.info("Reading configuration file")
    config = PluginConfig.from_xml(args.config_xml_file)
    plugin_elem = config.plugins[args.plugin_name]
    plugin_name = plugin_elem.get("name")
    plugin_func = PLUGINS[plugin_name]
    # read library table
    logging.info("Parsing library table")
    libs = []
    for lib in LibraryInfo.from_file(args.library_table):
        if not lib.is_valid():
            logging.warning("\t[SKIPPED] Library %s not valid" % (lib.library_id))
            continue
        libs.append(lib)
    logging.info("\tfound %d libraries" % (len(libs)))
    logging.info("Writing plugin scripts")
    retcode = plugin_func(libs, config, plugin_elem, 
                          args.keep_tmp, args.dryrun)
    logging.info("Done")
    return retcode

if __name__ == '__main__':
    sys.exit(main())