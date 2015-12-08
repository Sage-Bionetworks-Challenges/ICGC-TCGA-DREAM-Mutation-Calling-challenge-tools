#!/usr/bin/env python

import sys, os
import vcf

'''
Submission evaluation code for TCGA/ICGC/DREAM SMC
Adam Ewing, ewingad@soe.ucsc.edu
Requires PyVCF (https://github.com/jamescasbon/PyVCF)
'''

def match(subrec, trurec, vtype='SNV'):
    assert vtype in ('SNV', 'SV', 'INDEL')

    if vtype == 'SNV' and subrec.is_snp and trurec.is_snp:
        if subrec.POS == trurec.POS and subrec.REF == trurec.REF and subrec.ALT == trurec.ALT:
            return True

    if vtype == 'INDEL' and subrec.is_indel and trurec.is_indel:
        if subrec.POS == trurec.POS and subrec.REF == trurec.REF and subrec.ALT == trurec.ALT:
            return True

    if vtype == 'SV' and subrec.is_sv and trurec.is_sv:
        trustart, truend = expand_sv_ends(trurec, useCIs=False)
        substart, subend = expand_sv_ends(subrec, useCIs=False)

        # check for overlap
        if min(truend, subend) - max(trustart, substart) > 0:
            return True

    return False


def expand_sv_ends(rec, useCIs=True):
    '''
    assign start and end positions to SV calls
    using conf. intervals if present and useCIs=True
    '''
    startpos, endpos = rec.start, rec.end
    assert rec.is_sv

    try:
        if rec.INFO.get('END'): # sometimes this is a list, sometimes it's an int
            if isinstance(rec.INFO.get('END'), list):
                endpos = int(rec.INFO.get('END')[0])
            if isinstance(rec.INFO.get('END'), int):
                endpos = int(rec.INFO.get('END'))

        if useCIs:
            if rec.INFO.get('CIPOS'):
                ci = map(int, rec.INFO.get('CIPOS'))
                if ci[0] < 0:
                    startpos += ci[0]

            if rec.INFO.get('CIEND'):
                ci = map(int, rec.INFO.get('CIEND')) 
                if ci[0] > 0:
                    endpos += ci[0]

    except TypeError as e:
        sys.stderr.write("error expanding sv interval: " + str(e) + " for record: " + str(rec) + "\n")

    if startpos > endpos:
        endpos, startpos = startpos, endpos

    return startpos, endpos


def relevant(rec, vtype, ignorechroms):
    ''' Return true if a record matches the type of variant being investigated '''
    rel = (rec.is_snp and vtype == 'SNV') or (rec.is_sv and vtype == 'SV') or (rec.is_indel and vtype == 'INDEL')

    # 'ignore' types are always excluded
    if rec.INFO.get('SVTYPE'):
        rec_svtype_raw = rec.INFO.get('SVTYPE')
        if isinstance(rec_svtype_raw, list):
            # take the first value if SVTYPE is a list
            rec_svtype_raw = rec_svtype_raw[0]
        if rec_svtype_raw in ('IGN', 'MSK'):
            rel = False 

    return rel and (ignorechroms is None or rec.CHROM not in ignorechroms)


def passfilter(rec):
    ''' Return true if a record is unfiltered or has 'PASS' in the filter field (pyvcf sets FILTER to None) '''
    if rec.FILTER is None or rec.FILTER == '.' or not rec.FILTER:
        return True
    return False


def mask(rec, vcfh, truchroms, debug=False, active=True):
    ''' mask calls in IGN/MSK regions '''

    if rec.CHROM in truchroms:
        if rec.is_sv:
            teststart, testend = expand_sv_ends(rec, useCIs=False)
        else:
            teststart = rec.POS - 1
            testend = rec.POS

        for overlap_rec in vcfh.fetch(rec.CHROM, teststart, testend):
            overlap_rec_svtype = overlap_rec.INFO.get('SVTYPE')
            if overlap_rec_svtype and isinstance(overlap_rec_svtype, list):
                # take the first value if SVTYPE is a list
                overlap_rec_svtype = overlap_rec_svtype[0]
            
            # SNV and INDEL calls are ignored in IGN regions, cannot be deactivated with active=False 
            if (rec.is_snp or rec.is_indel) and overlap_rec_svtype == 'IGN':
                return True

            # calls are ignored in MSK regions if active=True
            if overlap_rec_svtype == 'MSK':
                is_masked = True
                
                if rec.is_sv:
                    mskstart, mskend = expand_sv_ends(overlap_rec, useCIs=False)
                    
                    # compute fraction of rec that is in the masked region
                    overlap_frac = float(min(mskend, testend) - max(mskstart, teststart))/float(testend - teststart)
                    if overlap_frac <= 0.5:
                        # sizable fraction of subrec is in unmasked region
                        is_masked = False
                
                if is_masked:
                    if debug:
                        print "DEBUG: submitted:", str(rec), "overlaps:", str(overlap_rec)
                    if active:
                        return True
                
    return False


def countrecs(result):
    ''' return number of counted mutation calls in submission '''
    assert 'tp' in result and 'fp' in result, "invalid result dictionary!"
    
    ncalls = result['tp'] + result['fp']
    
    return ncalls


def prefix(rec, usechr):
    ''' adjust presence/absence of "chr" prefix according to whether usechr is True or False '''
    if usechr and not rec.CHROM.startswith('chr'):
        rec.CHROM = 'chr' + rec.CHROM
    
    if not usechr and rec.CHROM.startswith('chr'):
        rec.CHROM = rec.CHROM.replace('chr', '')
    
    return rec
    
    
def evaluate(submission, truth, vtype='SNV', ignorechroms=None, truthmask=True):
    ''' return TP, FP and FN counts '''

    assert vtype in ('SNV', 'SV', 'INDEL')
    subvcfh = vcf.Reader(filename=submission)
    truvcfh = vcf.Reader(filename=truth)

    tpcount = 0 # counts TPs that are not masked
    fpcount = 0 # counts FPs that not masked
    tpcountmasked = 0 # counts masked TPs
    fpcountmasked = 0 # counts masked FPs
    subrecs = 0 # counts all predicted variant records
    subuncount = 0 # counts predicted variant records that are not counted as TPs/FPs
    trunotmasked = 0 # counts true variants that are not masked

    truchroms = {}

    ''' store list of truth records, otherwise the iterator needs to be reset '''
    trulist = [trurec for trurec in truvcfh]
    
    ''' track whether the truth uses the "chr" prefix (all truth entries are assumed to use the same reference) '''
    usechr = trulist[0].CHROM.startswith('chr')

    ''' count records in truth vcf, track contigs/chromosomes '''
    for trurec in trulist:
        if relevant(trurec, vtype, ignorechroms):
            truchroms[trurec.CHROM] = True
            if not mask(trurec, truvcfh, truchroms, active=truthmask):
                trunotmasked += 1
    
    # sanity check
    if trunotmasked == 0:
        raise Exception("No unmasked records found in truth file!\n")


    '''
    keep track of 'truth' sites used, they should only be usable once
    if toCount=True, the true variant will be counted; if False, it will not
    '''
    used_truth = {}
    
    '''
    if submitters use MATEID in their BND calls we can 'tie' them together,
    indexed by one mate, contains info on other mate in pair
    '''
    used_bnd_mates = {}

    ''' parse submission vcf, compare to truth '''
    for subrec in subvcfh:
        subrec = prefix(subrec, usechr)
        if relevant(subrec, vtype, ignorechroms) and passfilter(subrec):
            subrecs += 1
            
            matched = 'UN'

            startpos, endpos = subrec.start, subrec.end

            if vtype == 'SV' and subrec.is_sv:
                startpos, endpos = expand_sv_ends(subrec, useCIs=False)
            
            sub_is_masked = mask(subrec, truvcfh, truchroms, active=truthmask)
            
            if subrec.CHROM in truchroms:
                truoverlaplist = [trurec for trurec in truvcfh.fetch(subrec.CHROM, startpos, end=endpos)]
                for trurec in truoverlaplist:
                    if relevant(trurec, vtype, ignorechroms) and match(subrec, trurec, vtype=vtype):
                        # subrec matches a true variant
                        tru_is_masked = mask(trurec, truvcfh, truchroms, active=truthmask)
                        
                        if not sub_is_masked and not tru_is_masked:
                            if str(trurec) not in used_truth and matched == 'TP':
                                # subrec already matches another true variant
                                used_truth[str(trurec)] = { 'toCount': False, 'masked': tru_is_masked }
                            
                            elif matched != 'TP' and (str(trurec) not in used_truth or not used_truth[str(trurec)]['toCount']):
                                # subrec matches an unused true variant
                                matched = 'TP'
                                used_truth[str(trurec)] = { 'toCount': True, 'masked': tru_is_masked }
                        
                        elif str(trurec) not in used_truth:
                            # subrec matches an unused true variant but at least one is masked
                            used_truth[str(trurec)] = { 'toCount': False, 'masked': tru_is_masked }
            
                        if matched != 'TP':
                            '''
                            matched a true variant that was already counted
                            and/or has a different mask status/both are masked
                            only note if a true variant wasn't already IDed for subrecj
                            '''
                            matched = 'T'

            if matched == 'TP' or matched == 'T':
                if subrec.ID in used_bnd_mates and not used_bnd_mates[subrec.ID]['positive']:
                    # BND mate call was a false positive, remove conflict
                    if used_bnd_mates[subrec.ID]['masked']:
                        fpcountmasked -= 1
                    else:
                        fpcount -= 1
                    subuncount += 1
                
                elif subrec.INFO.get('MATEID'):
                    # keep track of the mate info
                    if isinstance(subrec.INFO.get('MATEID'), list):
                        mateID = subrec.INFO.get('MATEID')[0]
                    else:
                        mateID = subrec.INFO.get('MATEID')
                    used_bnd_mates[mateID] = { 'positive': True, 'masked': sub_is_masked }
            
            elif matched == 'UN' and subrec.ID not in used_bnd_mates:
                # subrec does not match a true variant and it is not tied to a previous match through a mate
                matched = 'FP'
                
                if sub_is_masked:
                    fpcountmasked += 1
                else:
                    fpcount += 1
                    
                if subrec.INFO.get('MATEID'):
                    # keep track of the mate info
                    if isinstance(subrec.INFO.get('MATEID'), list):
                        mateID = subrec.INFO.get('MATEID')[0]
                    else:
                        mateID = subrec.INFO.get('MATEID')
                    used_bnd_mates[mateID] = { 'positive': False, 'masked': sub_is_masked }
                
            if matched == 'TP':
                if sub_is_masked:
                    tpcountmasked += 1
                else:
                    tpcount += 1
            elif matched == 'UN' or matched == 'T':
                subuncount += 1

    # sanity checks
    assert (tpcount + fpcount + tpcountmasked + fpcountmasked + subuncount == subrecs)
    
    if subrecs == 0:
        raise Exception("No filter-passing variants in submission! Are you sure you selected the correct variant type (SNV/INDEL/SV)?\n")

    # count the true variants that should not be counted as TPs/FNs
    truuncountnotmasked = 0
    for k, v in used_truth.iteritems():
        if not v['toCount'] and not v['masked']:
            truuncountnotmasked += 1
    
    result = { 'tp' : float(tpcount),
               'fp' : float(fpcount),
               'fn' : float(trunotmasked - tpcount - truuncountnotmasked) }
    
    return result


def stats(result):
    ''' calculate precision, recall, fscore  from result dictionary '''
    assert 'tp' in result and 'fp' in result and 'fn' in result, "invalid result dictionary!"
    
    recall = float(1)
    if result['tp'] + result['fn'] > 0:
        recall = result['tp'] / (result['tp'] + result['fn'])
    
    precision = float(1)
    if result['tp'] + result['fp'] > 0:
        precision = result['tp'] / (result['tp'] + result['fp'])
    
    fscore = float(0)
    if precision > 0 or recall > 0:
        fscore = 2*((precision*recall) / (precision+recall))
    
    return recall, precision, fscore
    

if __name__ == '__main__':
    if len(sys.argv) == 4 or len(sys.argv) == 5:
        subvcf, truvcf, evtype = sys.argv[1:4]

        chromlist = None
        if len(sys.argv) == 5:
            chromlist = sys.argv[4].split(',')

        if not subvcf.endswith('.vcf') and not subvcf.endswith('.vcf.gz'):
            sys.stderr.write("submission VCF filename does not enc in .vcf or .vcf.gz\n")
            sys.exit(1)

        if not os.path.exists(truvcf + '.tbi'):
            sys.stderr.write("truth VCF does not appear to be indexed. bgzip + tabix index required.\n")
            sys.exit(1)

        if evtype not in ('SV', 'SNV', 'INDEL'):
            sys.stderr.write("last arg must be either SV, SNV, or INDEL\n")
            sys.exit(1)

        print "\nmasked:"
        counts = evaluate(subvcf, truvcf, vtype=evtype, ignorechroms=chromlist, truthmask=True)
        statresults = stats(counts)
        ncalls  = countrecs(counts)
        print "recall, precision, F1-score: " + ','.join(map(str, statresult))
        print "number of counted mutations in submission: " + str(ncalls)

        print "\nunmasked:"
        counts = evaluate(subvcf, truvcf, vtype=evtype, ignorechroms=chromlist, truthmask=False)
        statresults = stats(counts)
        ncalls  = countrecs(counts)
        print "recall, precision, F1-score: " + ','.join(map(str, statresult))
        print "number of counted mutations in submission: " + str(ncalls)

    else:
        print "standalone usage for testing:", sys.argv[0], "<submission VCF> <truth VCF (tabix-indexed)> <SV, SNV, or INDEL> [ignore chrom list (comma-delimited, optional)]"
