#! ${params.python3}

from Bio import SeqIO
from collections import defaultdict
import pickle
from ete3 import ncbi_taxonomy


def get_clade_dict(tax_map):
    ncbi = ncbi_taxonomy.NCBITaxa()
    taxon_clade_mod = defaultdict(str)
    for k, v in tax_map.items():
        lineage = ncbi.get_lineage(v)
        lineage_string = ""
        for l in lineage:
            rank = ncbi.get_rank([l])[l]
            if rank != 'no rank':
                lineage_string = "{} {};".format(lineage_string, ncbi.get_taxid_translator([l])[l])
        taxon_clade_mod[k] = lineage_string
    return taxon_clade_mod


def parse_eggnog_map(mapping_file):
    local_eggnog_map = {}
    with open(mapping_file) as handle:
        for line in handle:
            if line.startswith('-'):
                continue
            else:
                line = line.split()
                local_eggnog_map[line[1]] = line[0]
    return local_eggnog_map


def read_tax_map(fasta):
    tax_map = defaultdict(str)
    for rec in SeqIO.parse(fasta, 'fasta'):
        s = rec.id.split('.')
        tax_map[rec.id] = int(s[0])
    return tax_map


def commonprefix(taxon_paths):
    taxon_path_lists = [t.split('; ') for t in taxon_paths]
    s1 = min(taxon_path_lists)
    s2 = max(taxon_path_lists)
    for i, c in enumerate(s1):
        if c != s2[i]:
            return s1[:i]
    return s1


def clean_rec_id(recid):
    illegal_chars = ":,();[]' /+%\$@#&*<=>^!?\\{\\}\\"`"
    for c in illegal_chars:
        recid = recid.replace(c, '')
    return recid


def parse_fasta(fasta):
    seq_dict = defaultdict(list)
    seq_rec_dict = {}
    for rec in SeqIO.parse(fasta, 'fasta'):
        rec.id = clean_rec_id(rec.id)
        seq_dict[str(rec.seq)].append(rec.id)
        seq_rec_dict[rec.id] = rec
    return seq_rec_dict, seq_dict


def remove_duplicate_seqs(fasta, taxon_clade_mod, tax_map):
    ncbi = ncbi_taxonomy.NCBITaxa()
    seq_rec_dict, seq_dict = parse_fasta(fasta)
    removed_ids = {}
    for k,v in seq_dict.items():
        if len(v) > 1:
            tax_strings = []
            merged_ids = []
            seqrec = seq_rec_dict[v[0]]
            for id in v:
                tax_strings.append(taxon_clade_mod[id])
                merged_ids.append(id)
                taxon_clade_mod.pop(id)
                seq_rec_dict.pop(id)
                tax_map.pop(id)
            try:
                tax_path = commonprefix(tax_strings)
                taxon = tax_path[-1].replace(';', '')
                id = ncbi.get_name_translator([taxon])[taxon][0]
                seqrec.id = "{}.{}".format(id, clean_rec_id(taxon))
                removed_ids[seqrec.id] = merged_ids
                seqrec.description = ""
                taxon_clade_mod[seqrec.id] = "; ".join(tax_path)
                seq_rec_dict[seqrec.id] = seqrec
                tax_map[seqrec.id] = id
            except:
                print("key error")
    return taxon_clade_mod, seq_rec_dict, tax_map, removed_ids


def write_class_map(eggnog_id_mapping, basename):
    with open("{}.class".format(basename), 'w') as class_map:
        class_map.write("{}\\t{}\\n".format(eggnog_id_mapping, basename))


def write_unique_seqs(seq_rec_dict, basename):
    with open('{}.unique.fasta'.format(basename), 'w') as out:
        SeqIO.write(seq_rec_dict.values(), out, 'fasta')


def write_eggnog_map(eggnog_id_mapping, basename, tax_map):
    with open("{}.eggnog.map".format(basename), 'w') as eggnog_map:
        for k, v in tax_map.items():
            s = k.split('.')
            eggnog_map.write("{}\\t{}\\n".format(s[1], eggnog_id_mapping))


def write_taxon_map(taxon_clade_mod, basename):
    with open("{}.taxid.map".format(basename), 'w') as tax_map:
        for k,v in taxon_clade_mod.items():
            tax_map.write("{}\\t{}\\n".format(k,v))


def write_removed_seqs(removed_ids, basename):
    with open("{}.removed.seqs".format(basename), 'w') as removed_seq_out:
        for k,v in removed_ids.items():
            removed_seq_out.write("{}\\t{}\\n".format(k, "; ".join(v)))


def write_mapping_files(eggnog_id, local_eggnog_map, fasta):
    tax_map = read_tax_map(fasta)
    taxon_clade_mod = get_clade_dict(tax_map)
    taxon_clade_mod, seq_rec_dict, tax_map, removed_ids = remove_duplicate_seqs(fasta, taxon_clade_mod, tax_map)

    write_class_map(local_eggnog_map[eggnog_id], eggnog_id)
    write_unique_seqs(seq_rec_dict, eggnog_id)
    write_eggnog_map(local_eggnog_map[eggnog_id], eggnog_id, tax_map)
    write_taxon_map(taxon_clade_mod, eggnog_id)
    write_removed_seqs(removed_ids, eggnog_id)


def main(eggnog_id, eggnog_map, fasta):
    local_eggnog_map = parse_eggnog_map(eggnog_map)
    write_mapping_files(eggnog_id, local_eggnog_map, fasta)


if __name__ == '__main__':
    main("${fasta.simpleName}", "${megan_eggnog_map}", "${fasta}")
