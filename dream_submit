#!/usr/bin/env python

import sys
import os
import argparse
import json
import gzip
import re
import traceback
try:
    import synapseclient
    from synapseclient import File, Folder, Project
    from synapseclient import Evaluation, Submission, SubmissionStatus
except ImportError:
    print "Please Install Synapse Client Library"
    print ">>> pip install synapseclient"
    sys.exit(1)

try:
    import vcf
except ImportError:
    vcf = None

#Some of the evaluation interface methods require an up-to-date copy of the Synapse client
try:
    from distutils.version import StrictVersion
    if StrictVersion(re.sub(r'\.dev\d+', '', synapseclient.__version__)) < StrictVersion('0.5.1'):
        print "Please Upgrade Synapse Client Library"
        print ">>> pip install -U synapseclient"
        sys.exit(1)
except ImportError:
    pass


#The DOWNLOAD_ENTITY is a folder with several entities that represent tumor samples
#these entities are annotated with the evaluation IDs submissions should be sent to
DOWNLOAD_ENTITY="syn2280639"

# valid SVTYPE entries according to VCFv4.1
SV_TYPES = ('BND', 'CNV', 'DEL', 'DEL:ME', 'DUP', 'DUP:TANDEM', 'INS', 'INS: ME', 'INV')

'''
validatevcf

Attempts to validate a VCF that may contain SNVs, INDELs, and SVs
The main test is whether the VCF can be successfully parsed by PyVCF,
a few specific fields are check for as well (VCFs should define somatic,
imprecise breakends should have confidence intervals)

'''
def main_validate(args, syn):
    if vcf is None:
        print "Please instal PyVCF"
        print ">>> pip install PyVCF"
        sys.exit(1)
    validatevcf(args.vcf)

class ValidationError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def context(vcf, recnum):
    vcf_h = None
    ctbuf = 5
    if vcf.endswith('.gz'):
        vcf_h = gzip.open(vcf, 'rb')
    else:
        vcf_h = open(vcf, 'r')

    n = 0
    for line in vcf_h:
        if not line.startswith('#'):
            n += 1
            if n >= recnum-ctbuf and n <= recnum+ctbuf:
                print n,':',line.strip()
                found = True
            if n > recnum + ctbuf:
                vcf_h.close()
                return 
 
def validatevcf(invcf):
    assert invcf.endswith('.vcf') or invcf.endswith('.vcf.gz')
    print 'Starting Validation'
    try:
        vcfin = vcf.Reader(filename=invcf)
    except TypeError:
        if  invcf.endswith('.vcf.gz') and sys.version_info[:3] == (2,7,4):
            print "Error trying to open VCF file. There is a known bug in Python 2.7.4 that breaks PyVCF, please try a different version."
            sys.exit(0)
    recnum = 0
    try:
        indel_count = 0
        snv_count = 0
        sv_count = 0

        for rec in vcfin:
            ''' try to detemine how germline vs. somatic is specified '''
            somatic = False
            recnum += 1
            if rec.FILTER == 'GERMLINE' or rec.FILTER == 'SOMATIC':
                raise ValidationError('GERMLINE or SOMATIC in FILTER field: invalid VCF')
            if rec.INFO.get('SOMATIC'):
                somatic = True
            if str(rec.INFO.get('SS')).upper() == 'SOMATIC':
                somatic = True

            if rec.is_snp:
                snv_count += 1
            if rec.is_indel:
                indel_count += 1
            if rec.is_sv:
                sv_count += 1
                if rec.INFO.get('IMPRECISE'):
                    if not (rec.INFO.get('CIPOS')):
                        raise ValidationError('Imprecise SV record without CIPOS: ' + str(rec))

                    if rec.INFO.get('END') and not (rec.INFO.get('CIEND')):
                        raise ValidationError('Imprecise SV record using END without CIEND: ' + str(rec))
                
                if rec.INFO.get('SVTYPE'):
                    rec_svtype = rec.INFO.get('SVTYPE')
                    if isinstance(rec_svtype, list):
                        rec_svtype = rec_svtype[0]
                    
                    if rec_svtype not in SV_TYPES:
                        raise ValidationError('SV record with unrecognized SVTYPE: ' + rec_svtype)

            if not somatic:
                raise ValidationError('Record found not marked somatic, please mark somatic calls and do not include germline calls.')

        print "Validation complete."
        print "Total records:",recnum
        print "-"*60
        print "SNV count:",snv_count
        print "INDEL count:",indel_count
        print "SV count:",sv_count
        print "-"*60

    except:
        recnum += 1
        print "parse error in VCF on line",recnum
        print '-'*60
        traceback.print_exc(file=sys.stdout)
        print '-'*60
        print "context:"
        context(invcf, recnum)
        raise

def get_evaluation_by_name(syn, eval_name):
    if eval_name == 'test':
        ## Add test evaluations linked to the "MutationCallingChallengeStaging"
        ## project (syn2298928).
        out = {
            'snv':'2527266',
            'sv':'2527268',
            'indel':'2527270'
        }
        return out

    for field in ['uuid', 'sample_name']:
        out = None
        res = syn.query("select * from entity where parentId=='%s' and %s=='%s'" % (DOWNLOAD_ENTITY, field, eval_name) )
        for row in res['results']:
            if 'entity.snv_eval_id' in row or 'entity.sv_eval_id' in row or 'entity.indel_eval_id' in row:
                out = {}
                if 'entity.snv_eval_id' in row:
                    out['snv'] = row['entity.snv_eval_id'][0]
                if 'entity.sv_eval_id' in row:
                    out['sv'] = row['entity.sv_eval_id'][0]
                if 'entity.indel_eval_id' in row:
                    out['indel'] = row['entity.indel_eval_id'][0]
        if out is not None:
            return out
    return None


def get_evalutations(syn):
    out = {}
    res = syn.query("select * from entity where parentId=='%s' order by entity.name" % (DOWNLOAD_ENTITY) )
    for row in res['results']:
        if 'entity.uuid' in row:
            try:
                eval_map = {}
                if 'entity.snv_eval_id' in row:
                    eval_map['snv'] = row['entity.snv_eval_id'][0]
                if 'entity.sv_eval_id' in row:
                    eval_map['sv'] = row['entity.sv_eval_id'][0]
                if 'entity.indel_eval_id' in row:
                    eval_map['indel'] = row['entity.indel_eval_id'][0]
                out[row['entity.sample_name'][0]] = eval_map
            except KeyError:
                pass
    return out


def main_list_challenge(args, syn):
    e = get_evalutations(syn)
    print "SampleName\tsnv Evaluation\tsv Evaluation\tINDEL Evaluation"
    for sample_name, evaluation in e.items():
        if len(evaluation):
            print "%s\t%s\t%s\t%s" % (sample_name, evaluation.get('snv', "\t"), evaluation.get('sv', "\t"), evaluation.get('indel', "\t"))

def merge_lists(*lists):
    out = []
    for l in lists:
        for k in l:
            out.append(k)
    return out

def print_table(table):
    for row in table:
        print "|%s|" % " | ".join(row)

def main_list_sample(args, syn):
    res = syn.query("select * from entity where parentId=='%s'" % (DOWNLOAD_ENTITY))
    sample_map = {}
    for row in res['results']:
        if 'entity.sample_name' in row:
            #print  row['entity.sample_name'][0]
            file_name = row['entity.name']
            sample_name = row['entity.sample_name'][0]
            sample_type = row['entity.sample_type'][0]            
            uuid = row['entity.uuid'][0]
            file_type = row['entity.file_type'][0]
            file_md5 = row['entity.file_md5'][0]

            if sample_name not in sample_map:
                sample_map[sample_name] = {}
            if sample_type not in sample_map[sample_name]:
                sample_map[sample_name][sample_type] = {}
            if file_type not in sample_map[sample_name][sample_type]:
                sample_map[sample_name][sample_type][file_type] = []

            sample_map[sample_name][sample_type][file_type].append(
                (file_name, uuid, file_md5)
            )

            #print sample_name, sample_type, uuid, file_type
    
    for sample_name, sample_info in sample_map.items():
        if 'tumor' in sample_info and 'normal' in sample_info:
            print "##" + sample_name
            out = [ ['Tumor BAM Name', 'Tumor BAM UUID', 'Tumor BAM MD5'] ]
            out += sample_info['tumor']['bam']
            print_table(out)
            print ""
            if 'fastq' in  sample_info['tumor']:
                out = [ ['Tumor Fastq Name', 'Tumor Fastq UUID', 'Tumor Fastq MD5'] ]
                out += sample_info['tumor']['fastq']
                print_table(out)
                print ""

            out = [ ['Normal BAM Name', 'Normal BAM UUID', 'Normal BAM MD5'] ]
            out += sample_info['normal']['bam']
            print_table(out)
            print ""

            if 'fastq' in  sample_info['normal']:
                out = [ ['Normal Fastq Name', 'Normal Fastq UUID', 'Normal Fastq MD5'] ]
                out += sample_info['normal']['fastq']
                print_table(out)
                print ""
        

def main_join(args, syn):
    eval_ids = get_evaluation_by_name(syn, args.sample_name)
    if eval_ids is None:
        return

    if args.board == 'sv':
        evaluation = syn.getEvaluation(eval_ids['sv'])
    elif args.board == 'snv':
        evaluation = syn.getEvaluation(eval_ids['snv'])
    elif args.board == 'indel':
        evaluation = syn.getEvaluation(eval_ids['indel'])
    else:
        print "Error, please join 'sv' or 'snv' or 'indel', not " + args.board
        return
    try:
        print syn.joinEvaluation(evaluation)
    except synapseclient.exceptions.SynapseHTTPError:
        print "Error Trying to join Evaluation, did you already join?"


def main_submit(args, syn):
    """
    This method takes an existing VCF file, uploads it to a personal project folder, and then 
    submits it to the Dream Mutation calling Challenge
    """
    PUBLIC_SAMPLE_FOLDERS  = ["synthetic.challenge.set1.v2"]

    ## name sample folders nicely
    FOLDER_NAMES = {
        "synthetic.challenge.set1.v2": "IS1",
        "synthetic.challenge.set2":    "IS2",
        "synthetic.challenge.set3":    "IS3",
        "synthetic.challenge.set4":    "IS4",
        "synthetic.challenge.set5":    "IS5"
    }

    def folder_name(sample):
        return FOLDER_NAMES.get(sample, sample)

    def not_public(acl):
        new_permissions = [permissions for permissions in acl['resourceAccess'] \
            if 'principalId' in permissions \
            and permissions['principalId'] not in [synapseclient.PUBLIC, synapseclient.AUTHENTICATED_USERS]]
        acl['resourceAccess'] = new_permissions
        return acl


    if (args.vcf is None and args.entity is None) or args.name is None or args.team_name is None or args.sample_name is None or (args.vcf is not None and args.project_id is None):
        print """Usage:
dream_submit [sv|snv|indel] SAMPLE_NAME --vcf my_vcffile.gz --name "Name of Submission" --team-name "Team Name" --project-id syn12345
or
dream_submit [sv|snv|indel] SAMPLE_NAME --entity syn12345 --name "Name of Submission" --team-name "Team Name" 
"""
        if args.vcf is None:
            print "Please add --vcf or --entity"
        if args.sample_name is None:
            print "Please add SAMPLE_NAME"
        if args.name is None:
            print "Please add --name"
        if args.project_id is None:
            print "Please add --project-id"
        if args.team_name is None:
            print "Please add --team-name"
        sys.exit(0)
    
    eval_ids= get_evaluation_by_name(syn, args.sample_name)
    if eval_ids is None:
        print "ERROR: Not able to find evaluation. Check that the tumor sample name or uuid is correct."
        return

    if args.board == 'sv':
        evaluation = syn.getEvaluation(eval_ids['sv'])
    elif args.board == 'snv':
        evaluation = syn.getEvaluation(eval_ids['snv'])
    elif args.board == 'indel':
        evaluation = syn.getEvaluation(eval_ids['indel'])
    else:
        print "Error, please submit to 'sv' or 'snv' or 'indel', not " + args.board
        return

    if args.vcf is not None:
        ## require a valid VCF file for successful submission
        main_validate(args, None)

        ## place submission in a folder by sample name, that might already exist
        folder_id = syn._findEntityIdByNameAndParent(folder_name(args.sample_name), parent=args.project_id)
        if folder_id:
            ## if the folder's already there, don't mess with it
            folder = syn.get(folder_id)
        else:
            folder = Folder(folder_name(args.sample_name), parent=args.project_id, sample=args.sample_name)
            folder = syn.store(folder, createOrUpdate=True, forceVersion=False)
            ## If we create a folder, it's private unless it's in the
            ## set of public sample folders, in which case it inherits permission
            ## from the project. This preserves the participants' ability to keep
            ## things private, if they really want to.
            if args.sample_name not in PUBLIC_SAMPLE_FOLDERS:
                syn._storeACL(folder, not_public(syn._getACL(folder)))

        ## Subfolder by mutation type (snv, sv or indel), which is called "board" in the parameters
        ## folder structure will end up being Project/[sample]/[mutation_type]/[submission],
        ## for example: Project/IS1/snv/submission.txt
        subfolder = syn.store(Folder(args.board, parent=folder), createOrUpdate=True, forceVersion=False)

        entity = File(args.vcf, parent=subfolder)
        if args.test:
            entity.test_submission='true'
        syn.store(entity, forceVersion=False)
    else:
        entity = syn.get(args.entity)
        if args.test:
            entity.test_submission='true'
            entity = syn.store(entity, forceVersion=False)
    sub = syn.submit(evaluation, entity, name=args.name, teamName=args.team_name)

    print "You have submitted", sub


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Submit Files to the DREAM mutation calling challenge. Please see https://www.synapse.org/#!Synapse:syn312572/wiki/60703 for usage instructions.')
    #Stack.addJobTreeOptions(parser) 
    parser.add_argument("-u", "--user", help="UserName", default=None)
    parser.add_argument("-p", "--password", help="Password", default=None)

    subparsers = parser.add_subparsers(title="subcommand")

    parser_list_challenge = subparsers.add_parser('list-challenge')
    parser_list_challenge.set_defaults(login=True)
    parser_list_challenge.set_defaults(func=main_list_challenge)

    parser_list_sample = subparsers.add_parser('list-sample')
    parser_list_sample.set_defaults(login=True)
    parser_list_sample.set_defaults(func=main_list_sample)


    parser_join = subparsers.add_parser('join')
    parser_join.add_argument("board", choices=['snv', 'sv', 'indel'])
    parser_join.add_argument("sample_name", help="Sample name (or a UUID of the tumor file) this analysis was done on")    
    parser_join.set_defaults(login=True)
    parser_join.set_defaults(func=main_join)

    parser_submit = subparsers.add_parser('submit')
    parser_submit.add_argument("board", choices=['snv', 'sv', 'indel'])
    parser_submit.add_argument("sample_name", help="Sample name (or a UUID of the tumor file) this analysis was done on")    
    parser_submit.add_argument("--vcf", help="Path to the VCF file to be submitted")
    parser_submit.add_argument("--entity", help="Synapse Entity ID of VCF file to be submitted")
    parser_submit.add_argument("--name", required=True, help="Name of the submission")
    parser_submit.add_argument("--team-name", required=True, help="Name Team")
    parser_submit.add_argument("--project-id", default=None, help="The SYN id of your personal private working directory")
    parser_submit.add_argument("--test", action='store_true', default=False, help="Mark the submission as a test not to be included on the leaderboard")
    parser_submit.set_defaults(login=True)
    parser_submit.set_defaults(func=main_submit)

    indel_submit = subparsers.add_parser('indel')
    indel_submit.add_argument("sample_name", help="Sample name (or a UUID of the tumor file) this analysis was done on")    
    indel_submit.add_argument("--vcf", help="Path to the VCF file to be submitted")
    indel_submit.add_argument("--entity", help="Synapse Entity ID of VCF file to be submitted")
    indel_submit.add_argument("--name", required=True, help="Name of the submission")
    indel_submit.add_argument("--team-name", required=True, help="Name Team")
    indel_submit.add_argument("--project-id", default=None, help="The SYN id of your personal private working directory")
    indel_submit.add_argument("--test", action='store_true', default=False, help="Mark the submission as a test not to be included on the leaderboard")
    indel_submit.set_defaults(login=True)
    indel_submit.set_defaults(board="indel")
    indel_submit.set_defaults(func=main_submit)

    sv_submit = subparsers.add_parser('structural_variant')
    sv_submit.add_argument("sample_name", help="Sample name (or a UUID of the tumor file) this analysis was done on")    
    sv_submit.add_argument("--vcf", help="Path to the VCF file to be submitted")
    sv_submit.add_argument("--entity", help="Synapse Entity ID of VCF file to be submitted")
    sv_submit.add_argument("--name", required=True, help="Name of the submission")
    sv_submit.add_argument("--team-name", required=True, help="Name Team")
    sv_submit.add_argument("--project-id", default=None, help="The SYN id of your personal private working directory")
    sv_submit.add_argument("--test", action='store_true', default=False, help="Mark the submission as a test not to be included on the leaderboard")
    sv_submit.set_defaults(login=True)
    sv_submit.set_defaults(board="sv")
    sv_submit.set_defaults(func=main_submit)

    snv_submit = subparsers.add_parser('snv')
    snv_submit.add_argument("sample_name", help="Sample name (or a UUID of the tumor file) this analysis was done on")    
    snv_submit.add_argument("--vcf", help="Path to the VCF file to be submitted")
    snv_submit.add_argument("--entity", help="Synapse Entity ID of VCF file to be submitted")
    snv_submit.add_argument("--name", required=True, help="Name of the submission")
    snv_submit.add_argument("--team-name", required=True, help="Name Team")
    snv_submit.add_argument("--project-id", default=None, help="The SYN id of your personal private working directory")
    snv_submit.add_argument("--test", action='store_true', default=False, help="Mark the submission as a test not to be included on the leaderboard")
    snv_submit.set_defaults(login=True)
    snv_submit.set_defaults(board="snv")
    snv_submit.set_defaults(func=main_submit)

    parser_validate = subparsers.add_parser('validate')
    parser_validate.add_argument("vcf", help="VCF file to be validated")
    parser_validate.set_defaults(login=False) #don't need to log into Synapse to validate a file.
    parser_validate.set_defaults(func=main_validate)


    args = parser.parse_args() 
    if args.login:
        syn = synapseclient.Synapse()
        if args.user is not None:
            syn.login(args.user, args.password)
        else:
            if 'SYNAPSE_APIKEY' in os.environ and 'SYNAPSE_EMAIL' in os.environ:
                syn.login(email=os.environ['SYNAPSE_EMAIL'], apiKey=os.environ['SYNAPSE_APIKEY'])
            else:
                syn.login()
    else:
        syn = None

    sys.exit(args.func(args, syn))
