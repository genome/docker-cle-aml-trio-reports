FROM perl:5.20
MAINTAINER Feiyu Du <fdu@wustl.edu> 

LABEL \ 
  description="Perl scripts for AML Trio cwl reports" 

COPY full_variant_report.pl  /usr/local/bin/ 
COPY alignment_stat.pl /usr/local/bin/
COPY coverage_stat.pl /usr/local/bin/
