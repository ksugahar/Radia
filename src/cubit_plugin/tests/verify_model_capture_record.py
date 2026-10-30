"""Guards: actual C++ JSON decoding and source IDs. Assumes: synthetic native fixture."""
import json
import subprocess
import sys


def main():
    result = subprocess.run([sys.argv[1]], capture_output=True, check=True, timeout=30)
    data = json.loads(result.stdout.decode('utf8'))
    assert data['schema'] == 'cubit.model-capture.v1', data
    assert data['units'] == 'mm', data
    assert data['labels']['schema'] == 'cubit.source-face-labels.v1', data
    volume, = data['volumes']
    assert volume['id'] == 7 and volume['primitive']['center_mm'] == [1,2,3], volume
    assert volume['primitive']['radius_mm'] == 1, volume
    block, = data['labels']['blocks']
    assert block == dict(id=17,name='Body',source_volumes=[7],source_surface_ids=[11]), block
    sideset, = data['labels']['sidesets']
    assert sideset['name'] == 'Outer"wall' and sideset['source_blocks'] == [17], sideset
    assert sideset['source_surface_ids'] == [11], sideset
    assert data['physical_conditions_inferred'] is False, data
    assert data['geometry_accuracy_accepted'] is False, data
    print('Independent JSON decoder: source IDs, units, primitive and labels PASS')


if __name__ == '__main__':
    main()
