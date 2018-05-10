from lab.gene.anno import *
import pytest

GENIC_REGIONS = genic_region_list()


def test_get_genic_region_val():
    genic_region_to_bool = {genic_region: False for genic_region in GENIC_REGIONS}

    genic_region_to_bool['intergenic'] = True
    assert region_val_by_dict(genic_region_to_bool) == 0

    genic_region_to_bool['intergenic'] = False
    genic_region_to_bool['ORF'] = True
    genic_region_to_bool['intronic'] = True
    assert region_val_by_dict(genic_region_to_bool) == 258


def test_parse_genic_region_val():
    genic_region_to_bool = parse_genic_region_val(0)

    for genic_region in GENIC_REGIONS:
        if genic_region == 'intergenic':
            assert genic_region_to_bool[genic_region] is True
        else:
            assert genic_region_to_bool[genic_region] is False

    genic_region_to_bool = parse_genic_region_val(258)

    for genic_region in GENIC_REGIONS:
        if genic_region == 'ORF' or genic_region == 'intronic':
            assert genic_region_to_bool[genic_region] is True
        else:
            assert genic_region_to_bool[genic_region] is False


def test_parse_genic_region_val_fail():
    with pytest.raises(AssertionError):
        parse_genic_region_val(2 ** len(GENIC_REGIONS))

    with pytest.raises(AssertionError):
        parse_genic_region_val(-1)
