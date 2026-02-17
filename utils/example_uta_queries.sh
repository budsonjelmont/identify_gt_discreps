# Set search path for convenience
set search_path = "$user", public, uta_20210129b;

# Describe the genome-transcript alignments view 
\d+ tx_exon_aln_v

# Example queries to find genome-transcript discrepancies
select * from tx_exon_aln_v where tx_ac in ('NM_000059.3','NM_000090.3','NM_000257.2','NM_000384.2','NM_004329.2')
and cigar !~ '^[0-9]+=$' and alt_ac ~ 'NC.*';

select cigar, tx_ac, alt_ac, ord, (tx_end_i - tx_start_i) as tx_ex_len, (alt_end_i - alt_start_i) as alt_ex_len  
from tx_exon_aln_v where tx_ac = 'NM_001166053.1' and cigar !~ '^[0-9]+=$' and alt_ac = 'NC_000004.11' order by ord;

select concat(tx_ac,'|',tx_start_i,'|',tx_end_i,'|',alt_ac,'|',alt_start_i,'|',alt_end_i,'|',alt_strand,'|',alt_aln_method,'|',ord,'|',cigar,'|') as id , tx_ac as tx, alt_ac as chr from tx_exon_aln_v where cigar ~ '^[0-9]+=[0-9]' limit 50;
