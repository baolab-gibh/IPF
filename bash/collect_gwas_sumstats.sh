#!/bin/bash
# File: collect_gwas_sumstats.sh
# Author: Zhenhua Zhang
# E-mail: zhang_zhenhua@gibh.ac.cn
# Created: Jun 03, 2024
# Updated: Jun 03, 2024


ieu_gwas_baseurl=http://gwas-api.mrcieu.ac.uk/swagger.json
curl -X GET "http://gwas-api.mrcieu.ac.uk/status" -H  "accept: application/json" 
curl -X GET "http://gwas-api.mrcieu.ac.uk/gwasinfo" -H  "accept: application/json" 
