DATA_DIR=../data

mkdir -p $DATA_DIR

MANE_URL='https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.2/MANE.GRCh38.v1.2.summary.txt.gz'
MANE_OUTFILE_GZ=${DATA_DIR}/mane_grch38.txt.gz
MANE_OUTFILE=${MANE_OUTFILE_GZ%.*}
wget ${MANE_URL} -O ${MANE_OUTFILE_GZ}
gunzip -d ${MANE_OUTFILE_GZ}
