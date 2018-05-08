#!/usr/bin/env python3

"""Given a GRE ID, prints a list of all motif and their PSSMs
in JSON format
"""

import pymongo
import argparse
import json

MONGO_HOST = 'como'
MONGO_PORT = 27018
DESCRIPTION = """mtb_gre_motifs - retrieve motifs for GRE
"""

######################################################################
### Export
###############################

MEME_FILE_HEADER = """MEME version 4.12.0

ALPHABET= ACGT

strands: + -
"""

def write_pssm(outfile, element):
    """writes a single PSSM to the given file"""
    evalue = 0
    num_sites = 20
    motif_name = "%s-%d-%d" % (element['run'], element['cluster'], element['motif_num'])

    outfile.write('\nMOTIF %s\n' % motif_name)
    pssm_rows = element['pssm']

    outfile.write('letter-probability matrix: alength= 4 w= %d nsites= %d E= %.3e\n' % (len(pssm_rows), num_sites, evalue))
    for r in pssm_rows:
        outfile.write('%5.3f %5.3f %5.3f %5.3f\n' % (r['a'], r['c'], r['g'], r['t']))


if __name__ == '__main__':
    client = pymongo.MongoClient(host=MONGO_HOST, port=MONGO_PORT)
    db = client['mtu_db']
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('gre_id', type=int, help="GRE ID")
    parser.add_argument('outfile', help="path to output file")
    parser.add_argument('--format', help="output format (meme or json)", default="meme")
    args = parser.parse_args()

    result = db.motif_info.find({"gre_id": args.gre_id })
    if args.format == "json":
        data = []

    with open(args.outfile, 'w') as outfile:
        if args.format == 'meme':
            outfile.write(MEME_FILE_HEADER)
        for row in result:
            pssm = row['pwm']
            seqtype = row['seqtype']
            motif_num = row['motif_num']
            cluster_id = row['cluster_id']
            ri = db.bicluster_info.find({'_id': cluster_id})[0]
            run_id = ri['run_id']
            cluster_num = ri['cluster']
            run = db.ensemble_info.find({'run_id': run_id})[0]
            run_name = run['run_name']
            element = {
                "run": run_name,
                "cluster": cluster_num,
                "motif_num": motif_num,
                "seqtype": seqtype,
                "pssm": pssm
            }
            if args.format == 'json':
                data.append(element)
            else:
                write_pssm(outfile, element)
        if args.format == 'json':
            outfile.write(json.dumps(data))
