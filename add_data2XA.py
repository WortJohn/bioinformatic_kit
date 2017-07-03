#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division,print_function
import sys, re, getpass
import getpass, datetime
import yaml
import logging
import argparse
import database
import mymodule

def get_user_date():
    user_date = {}
    user_date['analyst'] = getpass.getuser()
    user_date['analysis_date'] = str(datetime.datetime.now())
    return user_date

def add_2db(all_infos):
    try:
        db2 = database.Db(db='database2')
        placeholders = ', '.join(['%s'] * len(all_infos))
        columns = ", ".join(all_infos.keys())
        insert_sql = "INSERT INTO XA (%s) VALUES (%s)" % (columns, placeholders)
        db2.cur.execute(insert_sql, all_infos.values())
        db2.conn.commit()
        db2.close()
    except Exception as e:
        logging.exception(e)

def cal_ot_dedup(qc_dict):
    qc_dict['ot'] = round(float(qc_dict['n_hit']) / float(qc_dict['n_map']), 3)
    qc_dict['dedup_rate'] = round(float(qc_dict['n_hit_dedup']) / float(qc_dict['n_hit']), 3)
    return qc_dict

def get_qc(qc):
    qc = yaml.load(open(qc))
    qc = cal_ot_dedup(qc)
    return qc

def read_probe_file(config, probe):
    if re.search(r'dna', probe, re.I):
        return config['DEFAULT']['dna_probe']
    else:
        return config['DEFAULT']['rna_probe']

def main(id, qc, probe, tab):
    try:
        probe_path = read_probe_file(mymodule.read_config_XA(), probe)
        cov_info = mymodule.cal_cov_percent(id, probe_path, tab)
        qc = get_qc(qc)
        all_info = dict(cov_info, **qc)
        all_info.update({'hybrid_probe':probe})
        all_info.update(get_user_date())
        add_2db(all_info)
    except Exception as e:
        logging.exception(e)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="The program is used for add QC info of blood cancer sample to MySQL database.")
    parser.add_argument('-i', dest='id', help='The analysis id. Required', required=True)
    parser.add_argument('-q', dest='qc_file', help='QC file. default: id.qc')
    parser.add_argument('-p', dest='probe_type', help='probe type which is used for hybrid capture', choices=['DNA', 'dna', 'RNA', 'rna'], required=True)
    parser.add_argument('-t', dest='tab_file', help='The snp tab file, such as: 2017.snp.vcf.tab. Default: id.snp.vcf.tab.')
    args = parser.parse_args()
    if not args.qc_file:
        args.qc_file = str(args.id) + ".qc"
    if not args.tab_file:
        args.tab_file = str(args.id) + ".snp.vcf.tab"
    main(args.id, args.qc_file, args.probe_type, args.tab_file)
