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

if __name__ == '__main__':
    client = pymongo.MongoClient(host=MONGO_HOST, port=MONGO_PORT)
    db = client['mtu_db']
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('gre_id', type=int, help="GRE ID")
    args = parser.parse_args()

    result = db.motif_info.find({"gre_id": args.gre_id })
    data = []
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
        data.append(element)
    print(json.dumps(data))
