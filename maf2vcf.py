import itertools
import vcf


def order(lst, ordering, key=None):
    """Sorts, in-place, and return lst sorted by ordering on key.
    Args:
        lst: a list to be sorted.
        ordering: a list defining an ordering on lst. All values/keyed values in
          lst must be in ordering.
    Optional Args:
        key: a string (name of attribute) or function which returns a value
          which will be ordered. If None, identity is used as the key.
    Use:
        order([132, 99, 22], ordering=[99, 22, 44, 132])
        # => [99, 22, 132]
        order([{'thing': 'bathrobe'}, {'thing': 'toga'}],
              ['toga', 'bathrobe'], key='thing')
        # => [{'thing': 'toga'}, {'thing': 'bathrobe'}]
    """
    if key is None:
        lookup = lambda x: x
    elif isinstance(key, basestring):
        lookup = lambda x: x[key]
    else: # the key is a function
        lookup = key
    ordering = {name: idx for idx, name in enumerate(ordering)}
    lst.sort(key=lambda x: ordering[lookup(x)])
    return lst

def vcf_format(val):
    """Format a value as stored in the database to a value as stored in a VCF.
    """
    if val is None:
        return '.'
    try:
        return int(val)
    except ValueError: pass
    try:
        return float(val)
    except ValueError: pass
    vals = val.split(',')
    if len(vals) <= 1:
        return val
    elif val[0] != '[':
        return ','.join(v.strip() for v in vals)
    else: # It's a list.
        assert val[0] == '['
        assert val[-1] == ']'
        return ','.join(v.strip() for v in vals)[1:-1]

def _call_key(gt):
    """Return a string key on the CHROM/POS/REF/ALT used for grouping."""
    return '{}-{}-{}-{}'.format(gt['contig'], gt['position'],
                                gt['reference'], gt['alternates'])


def genotypes_to_records(genotypes, reader, extant_columns):
    """Return a list of vcf.model._Record, converted from a list of genotype
    relations, taken from the database.
    This is used to write to a .vcf file with vcf.parser.Writer.
    Args:
        genotypes: a list of genotype relations from the database.
        reader: a template vcf.Reader.
        extant_columns: a list of columns which have values for these genotypes.
    """
    info_fields, format_fields = _fields_from_columns(extant_columns)
    CallData = vcf.model.make_calldata_tuple(format_fields)
    grouped_genotypes = itertools.groupby(genotypes, _call_key)
    sample_ordering = reader.samples
    records = []
    for _, gts in grouped_genotypes:
        samples = []
        gts = list(gts)
        gts = order(gts, sample_ordering, 'sample_name')
        for gt in gts:
            data = CallData(*[vcf_format(gt['sample:' + d]) for d in format_fields])
            call = vcf.model._Call(None, gt['sample_name'], data)
            samples.append(call)
        record = _make_record_from_gt(gts[0], info_fields, format_fields, samples)
        records.append(record)
    return records


def _fields_from_columns(columns):
    # 7 == 'sample:', 5 == 'info:' -- we're stripping them off the column names.
    format_fields = [c[7:] for c in columns if c.startswith('sample:')]
    info_fields = [c[5:] for c in columns if c.startswith('info:')]
    return info_fields, format_fields


def _maybe_split(string, char):
    if string:
        return string.split(char)
    else:
        return None


def _make_record_from_gt(genotype, info_fields, format_fields, samples):
    gt = genotype
    info = {name: vcf_format(gt['info:' + name]) for name in info_fields}
    return vcf.model._Record(gt['contig'],                         # CHROM
                             gt['position'],                       # POS
                             gt['id'],                             # ID
                             gt['reference'],                      # REF
                             _maybe_split(gt['alternates'], ','),  # ALT
                             gt['quality'],                        # QUAL
                             _maybe_split(gt['filters'], ','),     # FILTER
                             info,                                 # INFO
                             ':'.join(format_fields),              # FORMAT
                             [0],                               # sample_indexes
                             samples=samples)                   # samples/calls


def genotypes_to_file(genotypes, header, extant_columns, fd):
    """Write genotypes to a VCF file.
    Args:
        genotypes: a list of genotype relations from the database.
        extant_columns: a list of columns which have values for these genotypes.
        header: the text header of the original VCF file.
        fd: an open file which will be written to.
    """
    template = vcf.Reader(l for l in header.split('\n'))
    records = genotypes_to_records(genotypes, template, extant_columns)
    writer = vcf.Writer(fd, template)
    for record in records:
        writer.write_record(record)
    fd.seek(0, 0)

if __name__ == '__main__':
    header = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR'
    template = vcf.Reader(header)
    genotypes = \
    [
        {
            'sample_name': 'TUMOR',
            'contig': 'chr7',
            'position': 151879673,
            'id': '',
            'reference': 'G',
            'alternates': 'A',
            'quality': '',
            'filters': 'PASS',
            'sample:GT': '0/1'
        },
        {
            'sample_name': 'NORMAL',
            'contig': 'chr7',
            'position': 151879673,
            'id': '',
            'reference': 'G',
            'alternates': 'A',
            'quality': '',
            'filters': 'PASS',
            'sample:GT': '0/0'
        }
    ]
    extant_columns = ['sample:GT']
    outfile = open('/Users/hammer/Desktop/blah.vcf', 'w')
    genotypes_to_file(genotypes, header, extant_columns, outfile)


