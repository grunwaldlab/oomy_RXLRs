#!/usr/bin/env python
"""Implements assorted RXLR motif methods from the literature.

# Original author: Peter J. Cock
# Modified by: Nicholas C. Cauldron

This script takes exactly four command line arguments:
 * Protein FASTA filename
 * Number of threads
 * Model name (Bhattacharjee2006, Win2007, Whisson2007, effectorp3, secretome)
 * Output tabular filename

The model names are:
 * Bhattacharjee2006: Simple regular expression search for RXLR
   with additional requirements for positioning and signal peptide.
 * Win2007: Simple regular expression search for RXLR, but with
   different positional requirements.
 * Whisson2007: As Bhattacharjee2006 but with a more complex regular
   expression to look for RXLR-EER domain, and additionally calls HMMER.
 * effectorp3: Does not use regex, it runs EffectorP v3.
   Implemented in Cauldron et al. 2023
 * Secretome: Just runs SignalP on all candidates. Does not test any RXLR models.

See the help text in the accompanying Galaxy tool XML file for more
details including the full references.

Note
----
Bhattacharjee et al. (2006) and Win et al. (2007) used SignalP v2.0,
which is no longer available. The current release is SignalP v3.0
(Mar 5, 2007). We have therefore opted to use the NN Ymax position for
the predicted cleavage site, as this is expected to be more accurate.
Also note that the HMM score values have changed from v2.0 to v3.0.
Whisson et al. (2007) used SignalP v3.0 anyway.

Whisson et al. (2007) used HMMER 2.3.2, and althought their HMM model
can still be used with hmmsearch from HMMER 3, sadly this does give
slightly different results. We expect the hmmsearch from HMMER 2.3.2
(the last stable release of HMMER 2) to be present on the path under
the name hmmsearch2 (allowing it to co-exist with HMMER 3).

If using Conda, you should therefore install the special "hmmer2"
package from BioConda which provides "hmmsearch2" etc::

    conda install -c bioconda hmmer2

See https://bioconda.github.io/recipes/hmmer2/README.html and
https://anaconda.org/bioconda/hmmer2
"""

from __future__ import print_function

import os
import re
import subprocess
import sys
from time import gmtime, strftime

from seq_analysis_utils import fasta_iterator

if "-v" in sys.argv:
    print("RXLR Motifs v0.0.16")
    sys.exit(0)

if len(sys.argv) != 5:
    sys.exit(
        "Requires four arguments: protein FASTA filename, threads, "
        "model, and output filename"
    )

fasta_file, threads, model, tabular_file = sys.argv[1:]
hmm_output_file = tabular_file + ".hmm.tmp"
effp_output_file = tabular_file + ".effp.tsv"
effp_stdout_file = tabular_file + ".effp.tmp"
signalp_input_file = tabular_file + ".fasta.tmp"
signalp_output_file = tabular_file + "sp3_tabular.tmp"
min_signalp_hmm = 0.9
hmmer_search = "hmmsearch2"
effp = "~/opt/EffectorP-3.0/EffectorP.py"

if model.lower() == "bhattacharjee2006":
    signalp_trunc = 70
    re_rxlr = re.compile("R.LR")
    min_sp = 10
    max_sp = 40
    max_sp_rxlr = 100
    min_rxlr_start = 1
    # Allow signal peptide to be at most 40aa, and want RXLR to be
    # within 100aa, therefore for the prescreen the max start is 140:
    max_rxlr_start = max_sp + max_sp_rxlr
elif model.lower() == "win2007":
    signalp_trunc = 70
    re_rxlr = re.compile("R.LR")
    min_sp = 10
    max_sp = 40
    min_rxlr_start = 30
    max_rxlr_start = 60
    # No explicit limit on separation of signal peptide clevage
    # and RXLR, but shortest signal peptide is 10, and furthest
    # away RXLR is 60, so effectively limit is 50.
    max_sp_rxlr = max_rxlr_start - min_sp + 1
elif model.lower() == "whisson2007":
    signalp_trunc = 0  # zero for no truncation
    re_rxlr = re.compile("R.LR.{,40}[ED][ED][KR]")
    min_sp = 10
    max_sp = 40
    max_sp_rxlr = 100
    min_rxlr_start = 1
    max_rxlr_start = max_sp + max_sp_rxlr
elif model.lower() in {"effectorp3", "secretome"}:
    signalp_trunc = 70
    #re_rxlr = re.compile("R.LR")
    # Don't req any motif - can see overlap w/Bhatta for this
    re_rxlr = re.compile(".")
    min_sp = 10
    max_sp = 40
    max_sp_rxlr = 100
    min_rxlr_start = 1
    # Allow signal peptide to be at most 40aa, and want RXLR to be
    # within 100aa, therefore for the prescreen the max start is 140:
    max_rxlr_start = max_sp + max_sp_rxlr
else:
    sys.exit(
        "Did not recognise the model name %r\n"
        "Use Bhattacharjee2006, Win2007, Whisson2007, effectorp3, or secretome" % model
    )

print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))

def get_hmmer_version(exe, required=None):
    try:
        child = subprocess.Popen(
            [exe, "-h"],
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except OSError:
        raise ValueError("Could not run %s" % exe)
    stdout, stderr = child.communicate()
    if required:
        return required in stdout
    elif "HMMER 2" in stdout:
        return 2
    elif "HMMER 3" in stdout:
        return 3
    else:
        raise ValueError("Could not determine version of %s" % exe)

def get_effp_version(exe, required=None):
    try:
        child = subprocess.Popen(
            [exe, "-h"],
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except OSError:
        raise ValueError("Could not run %s" % exe) 
    stdout, stderr = child.communicate()

# Run hmmsearch for Whisson et al. (2007)
if model.lower() == "whisson2007":
    hmm_file = os.path.join(
        os.path.split(sys.argv[0])[0], "whisson_et_al_rxlr_eer_cropped.hmm"
    )
    if not os.path.isfile(hmm_file):
        sys.exit("Missing HMM file for Whisson et al. (2007)")
    if not get_hmmer_version(hmmer_search, "HMMER 2.3.2 (Oct 2003)"):
        sys.exit("Missing HMMER 2.3.2 (Oct 2003) binary, %s" % hmmer_search)

    hmm_hits = set()
    valid_ids = set()
    for title, seq in fasta_iterator(fasta_file):
        name = title.split(None, 1)[0]
        if name in valid_ids:
            sys.exit("Duplicated identifier %r" % name)
        else:
            valid_ids.add(name)
    if not valid_ids:
        # Special case, don't need to run HMMER if there are no sequences
        pass
    else:
        # I've left the code to handle HMMER 3 in situ, in case
        # we revisit the choice to insist on HMMER 2.
        hmmer3 = 3 == get_hmmer_version(hmmer_search)
        # Using zero (or 5.6?) for bitscore threshold
        if hmmer3:
            # The HMMER3 table output is easy to parse
            # In HMMER3 can't use both -T and -E
            cmd = "%s -T 0 --tblout %s --noali %s %s > /dev/null" % (
                hmmer_search,
                hmm_output_file,
                hmm_file,
                fasta_file,
            )
        else:
            # For HMMER2 we are stuck with parsing stdout
            # Put 1e6 to effectively have no expectation threshold (otherwise
            # HMMER defaults to 10 and the calculated e-value depends on the
            # input FASTA file, and we can loose hits of interest).
            cmd = "%s -T 0 -E 1e6 %s %s > %s" % (
                hmmer_search,
                hmm_file,
                fasta_file,
                hmm_output_file,
            )
        print("Running hmmsearch")
        return_code = os.system(cmd)
        if return_code:
            sys.stderr.write("Error %i from hmmsearch:\n%s\n" % (return_code, cmd))
            sys.exit(return_code)

        handle = open(hmm_output_file)
        for line in handle:
            if not line.strip():
                # We expect blank lines in the HMMER2 stdout
                continue
            elif line.startswith("#"):
                # Header
                continue
            else:
                name = line.split(None, 1)[0]
                # Should be a sequence name in the HMMER3 table output.
                # Could be anything in the HMMER2 stdout.
                if name in valid_ids:
                    hmm_hits.add(name)
                elif hmmer3:
                    sys.exit("Unexpected identifer %r in hmmsearch output" % name)
        handle.close()
        # if hmmer3:
        #     print "HMMER3 hits for %i/%i" % (len(hmm_hits), len(valid_ids))
        # else:
        #     print "HMMER2 hits for %i/%i" % (len(hmm_hits), len(valid_ids))
        # print "%i/%i matched HMM" % (len(hmm_hits), len(valid_ids))
        os.remove(hmm_output_file)
    del valid_ids

# Run effectorp3 if it's the selected model
if model.lower() == "effectorp3":
    # Check if we can run EffectorP3
    #if not get_effp_version(effp):
    #    sys.exit("Missing EffectorP 3 binary, %s" % effp)
    cyto_hits = set()
    apo_hits = set()
    non_hits = set()
    valid_ids = set()
    for title, seq in fasta_iterator(fasta_file):
        name = title.split(None, 1)[0]
        #name = title.split("\t")[0]
        if name in valid_ids:
            sys.exit("Duplicated identifier %r" % name)
        else:
            valid_ids.add(name)
    if not valid_ids:
        # Special case, don't need to run effectorp3 if there are no sequences
        pass
    else:
        cmd = "python3 %s -i %s -o %s > %s" % (
            effp,
            fasta_file,
            effp_output_file,
            effp_stdout_file
        )
        print("Running EffectorP3")
        print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        # Only run if there's not a file yet
        if not os.path.exists(effp_output_file):
            print("EffectorP output not found, running tool.")
            return_code = os.system(cmd)
            if return_code:
                sys.stderr.write("Error %i from EffectorP3:\n%s\n" % (return_code, cmd))
                sys.exit(return_code)

        handle = open(effp_output_file)
        for line in handle:
            #if not line.strip():
                # We do not expect blank lines in output table
            #    continue
            #elif line.startswith("#"):
            if line.startswith("#"):
                # Header
                continue
            else:
                hit = line.split("\t")
                name = line.split(None, 1)[0]
                cyto = hit[1]
                apo = hit[2]
                non_eff = hit[3]
                # Should be a sequence name in the HMMER3 table output.
                # Could be anything in the HMMER2 stdout.
                #print(name)
                #print(valid_ids)
                if name in valid_ids:
                    #hmm_hits.add(name)
                    # Check if classified as effector
                    if "Y" in cyto:
                        #print("Cyto hit %s: %s" % (name, cyto))
                        cyto_hits.add(name)
                    if "Y" in apo:
                        #print("Apo hit %s: %s" % (name, apo))
                        apo_hits.add(name)
                    if "Y" in non_eff:
                        #print("Non-eff hit %s: %s" % (name, non_eff))
                        non_hits.add(name)
                #elif hmmer3:
                else:
                    sys.exit("Unexpected identifer %r in EffectorP3 output" % name)
        handle.close()
        # if hmmer3:
        #     print "HMMER3 hits for %i/%i" % (len(hmm_hits), len(valid_ids))
        # else:
        #     print "HMMER2 hits for %i/%i" % (len(hmm_hits), len(valid_ids))
        # print "%i/%i matched HMM" % (len(hmm_hits), len(valid_ids))
        try:
            os.remove(effp_stdout_file)
        except:
            pass
    del valid_ids
    #print("Cytoplasmic effectors:")
    #print(cyto_hits)
    #print("Apoplastic effectors:")
    #print(apo_hits)

# Prepare short list of candidates containing RXLR to pass to SignalP
assert min_rxlr_start > 0, "Min value one, since zero based counting"
count = 0
total = 0
handle = open(signalp_input_file, "w")
for title, seq in fasta_iterator(fasta_file):
    total += 1
    name = title.split(None, 1)[0]
    # can't search regex with no criterion
    #if model != "effectorp3":
    #    match = re_rxlr.search(seq[min_rxlr_start - 1 :].upper())
    #elif model == "effectorp3":
    #    match = seq
    # Just keep original method b/c the string I set for effectorP hits all
    match = re_rxlr.search(seq[min_rxlr_start - 1 :].upper())
    # effectorp3 predicts cytoplasmic and apoplastic effectors so don't use regex
    # run effectorp3 for oomycetes
    if match and min_rxlr_start - 1 + match.start() + 1 <= max_rxlr_start:
        # This is a potential RXLR, depending on the SignalP results.
        # Might as well truncate the sequence now, makes the temp file smaller
        if signalp_trunc:
            handle.write(">%s (truncated)\n%s\n" % (name, seq[:signalp_trunc]))
        else:
            # Does it matter we don't line wrap?
            handle.write(">%s\n%s\n" % (name, seq))
        count += 1
handle.close()
print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
print("Running SignalP on %i/%i potentials." % (count, total))


# Run SignalP (using our wrapper script to get multi-core support etc)
signalp_script = os.path.join(os.path.split(sys.argv[0])[0], "signalp3.py")
if not os.path.isfile(signalp_script):
    sys.exit("Error - missing signalp3.py script")
cmd = "python '%s' 'euk' '%i' '%s' '%s' '%s'" % (
    signalp_script,
    signalp_trunc,
    threads,
    signalp_input_file,
    signalp_output_file,
)
return_code = os.system(cmd)
if return_code:
    sys.exit("Error %i from SignalP:\n%s" % (return_code, cmd))
print("SignalP done")
# potentially finish here if no model specified
if model.lower() == "secretome":
    print("Running SignalP without searching for effector motifs.")
    os.remove(signalp_input_file)
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    signalp_results = parse_signalp(signalp_output_file)
    total = 0
    sp_count_hmm = 0
    sp_count_nn_len = 0
    no_sp_count = 0
    for sp_hits in signalp_results:
        sp_id = sp_hits[0]
        sp_hmm_score = sp_hits[1]
        sp_nn_len = sp_hits[2]
        total += 1
        if sp_hmm_score >= min_signalp_hmm:
            sp_count_hmm += 1
            if min_sp <= sp_nn_len <= max_sp:
                sp_count_nn_len += 1
        else:
            no_sp_count += 1
    print("Out of %d proteins, %d passed HMM and %d had signal peptides between %d and %d" % (total, sp_count_hmm, sp_count_nn_len, min_sp, max_sp))
    print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    sys.exit()

def parse_signalp(filename):
    """Parse SignalP output, yield tuples of values.

    Returns tuples of ID, HMM_Sprob_score and NN predicted signal
    peptide length.

    For signal peptide length we use NN_Ymax_pos (minus one).
    """
    handle = open(filename)
    line = handle.readline()
    assert line.startswith("#ID\t"), line
    for line in handle:
        parts = line.rstrip("\t").split("\t")
        assert len(parts) == 20, repr(line)
        yield parts[0], float(parts[18]), int(parts[5]) - 1
    handle.close()

# Parse SignalP results and apply the strict RXLR criteria
total = 0
tally = {}
handle = open(tabular_file, "w")
handle.write("#ID\t%s\n" % model)
print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
print("Parse SignalP results and check for regex to assign final Y/N")
signalp_results = parse_signalp(signalp_output_file)
for title, seq in fasta_iterator(fasta_file):
    total += 1
    # Default is Non-effector
    rxlr = "N"
    name = title.split(None, 1)[0]
    match = re_rxlr.search(seq[min_rxlr_start - 1 :].upper())
    if match and min_rxlr_start - 1 + match.start() + 1 <= max_rxlr_start:
        del match
        # This was the criteria for calling SignalP,
        # so it will be in the SignalP results.
        sp_id, sp_hmm_score, sp_nn_len = next(signalp_results)
        assert name == sp_id, "%s vs %s" % (name, sp_id)
        # Only check regex for secreted proteins
        if sp_hmm_score >= min_signalp_hmm and min_sp <= sp_nn_len <= max_sp:
            match = re_rxlr.search(seq[sp_nn_len:].upper())
            # Check if secreted prot is also a regex match - if so, then convert to Y
            if match and match.start() + 1 <= max_sp_rxlr:  # 1-based counting
                rxlr_start = sp_nn_len + match.start() + 1
                if min_rxlr_start <= rxlr_start <= max_rxlr_start:
                    rxlr = "Y"
                    # For effp3 this just means passes signalp,
                    # so now check if protein passed our model
                    if model == "effectorp3":
                        #print("Checking effectorp 3 results: %s" % (name))
                        # Four way classifier: Y, both, apo, or neither
                        # count will be number of Y (cyto only), but include all the data
                        if name in apo_hits:
                            #print("apoplastic effector found")
                            rxlr = "apo"
                        if name in cyto_hits:
                            #print("cytoplasmic effector found")
                            #if name in apo_hits:
                            if rxlr == "apo":
                                #print("effector both cyto and apoplastic")
                                rxlr = "both"
                            else:
                                rxlr = "Y"
                        if name in non_hits:
                            #print("Non-effector")
                            rxlr = "N"
        #else:
        #    print("%s not secreted %s" % (name, sp_hmm_score))
    if model == "Whisson2007":
        # Combine the signalp with regular expression heuristic and the HMM
        if name in hmm_hits and rxlr == "N":
            rxlr = "hmm"  # HMM only
        elif rxlr == "N":
            rxlr = "neither"  # Don't use N (no)
        elif name not in hmm_hits and rxlr == "Y":
            rxlr = "re"  # Heuristic only
        # Now have a four way classifier: Y, hmm, re, neither
        # and count is the number of Y results (both HMM and heuristic)
    # Does this overwrite signalp results? Would explain why data is huge
    handle.write("%s\t%s\n" % (name, rxlr))
    try:
        tally[rxlr] += 1
    except KeyError:
        tally[rxlr] = 1
handle.close()
assert sum(tally.values()) == total

# Check the iterator is finished
try:
    next(signalp_results)
    assert False, "Unexpected data in SignalP output"
except StopIteration:
    pass

# Cleanup
os.remove(signalp_input_file)
#os.remove(signalp_output_file)

# Short summary to stdout for Galaxy's info display
print("%s for %i sequences:" % (model, total))
print(", ".join("%s = %i" % kv for kv in sorted(tally.items())))
print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
