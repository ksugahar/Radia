"""Guards: independent native/BDF/skin binding and explicit linear-mesh units.

Assumes: synthetic plane-defined skin; no licensed acquisition or app connection.
"""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SOURCE = Path(__file__).resolve().parents[1] / 'src/cubit_mesh_export/mesh_bundle.py'
spec = importlib.util.spec_from_file_location('pure_mesh_bundle_contract', SOURCE)
bundle = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bundle)


def fixture():
    # Outward skin from the six planes of an independently authored box.
    xyz = [(0.,0.,0.), (2.,0.,0.), (2.,3.,0.), (0.,3.,0.),
           (0.,0.,5.), (2.,0.,5.), (2.,3.,5.), (0.,3.,5.)]
    faces = [(1,4,3,2), (5,6,7,8), (1,2,6,5),
             (2,3,7,6), (3,4,8,7), (4,1,5,8)]
    return {'schema': 'cubit.native-mesh-capture.v1', 'units': 'mm',
            'nodes': [{'id': i, 'xyz': list(p)} for i,p in enumerate(xyz,1)],
            'elements': [{'source_id': 71, 'type': 'HEX8', 'region_id': 19,
                          'node_ids': list(range(1,9))}],
            'regions': [{'id': 19, 'name': 'box_region'}],
            'faces': [{'source_surface_id': i, 'source_face_id': i,
                       'source_volume_id': 2, 'source_volume_sense': 'Forward',
                       'node_ids': list(nodes)} for i,nodes in enumerate(faces,1)]}


def native_file(path, capture):
    # Independent standard NetCDF writer; never the production decoder.
    from scipy.io import netcdf_file
    with netcdf_file(path, 'w', version=2) as db:
        for name, size in [('num_nodes',8), ('num_elem',1), ('num_el_blk',1),
                           ('num_el_in_blk1',1), ('num_nod_per_el1',8)]:
            db.createDimension(name,size)
        db.createVariable('node_num_map','i',('num_nodes',))[:] = list(range(1,9))
        db.createVariable('elem_num_map','i',('num_elem',))[:] = [71]
        db.createVariable('eb_prop1','i',('num_el_blk',))[:] = [19]
        for i,axis in enumerate('xyz'):
            db.createVariable('coord'+axis,'d',('num_nodes',))[:] = [n['xyz'][i] for n in capture['nodes']]
        v = db.createVariable('connect1','i',('num_el_in_blk1','num_nod_per_el1'))
        v.elem_type = b'HEX8'
        v[:] = [list(range(1,9))]


class MeshBundleContracts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='mesh-bundle-contract-')
        self.addCleanup(self.tmp.cleanup)
        self.native = Path(self.tmp.name)/'alternate.e'
        self.mesh = fixture()
        native_file(self.native,self.mesh)

    def test_native_and_bdf_reparse_preserve_order_and_labels(self):
        nodes,cells = bundle.check_capture(self.mesh,self.native)
        data,mapping = bundle._write_bdf(nodes,cells,[19])
        self.assertIn(b'CHEXA,1,19,1,2,3,4,5,6,+N1',data)
        self.assertEqual(mapping,[{'source_id':71,'source_type':'HEX8','eid':1,'pid':19}])
        self.assertEqual(self.mesh['regions'][0]['name'],'box_region')

    def test_native_connectivity_mismatch_is_rejected(self):
        from scipy.io import netcdf_file
        with netcdf_file(self.native,'a') as db:
            db.variables['connect1'][:] = [[2,3,4,1,6,7,8,5]]
        with self.assertRaisesRegex(ValueError,'Independent native'):
            bundle.check_capture(self.mesh,self.native)

    def test_oriented_skin_cannot_be_self_reconstructed(self):
        for mode in ('reverse','missing'):
            mesh=copy.deepcopy(self.mesh)
            if mode=='reverse':mesh['faces'][0]['node_ids'].reverse()
            else:mesh['faces'].pop()
            with self.subTest(mode=mode),self.assertRaisesRegex(ValueError,'outward skin'):
                bundle.check_capture(mesh,self.native)

    def test_explicit_units_no_order_reduction_or_duplicate_ids(self):
        for field,value in [('units',''),('schema','unknown')]:
            mesh=copy.deepcopy(self.mesh);mesh[field]=value
            with self.assertRaises(ValueError):bundle.check_capture(mesh,self.native)
        mesh=copy.deepcopy(self.mesh);mesh['elements'][0]['node_ids'].append(9)
        with self.assertRaisesRegex(ValueError,'connectivity'):bundle.check_capture(mesh,self.native)
        mesh=copy.deepcopy(self.mesh);mesh['nodes'].append(mesh['nodes'][0])
        with self.assertRaisesRegex(ValueError,'node'):bundle.check_capture(mesh,self.native)

    def test_invalid_native_index_does_not_wrap_python_list(self):
        from scipy.io import netcdf_file
        with netcdf_file(self.native,'a') as db:db.variables['connect1'][0,0]=0
        with self.assertRaisesRegex(ValueError,'index outside'):bundle.check_capture(self.mesh,self.native)

    def test_twisted_hex_and_unknown_skin_sense_fail(self):
        mesh=copy.deepcopy(self.mesh);mesh['nodes'][6]['xyz']=[-2.,-3.,-5.]
        with self.assertRaisesRegex(ValueError,'Jacobian'):bundle.check_capture(mesh,self.native)
        mesh=copy.deepcopy(self.mesh);mesh['faces'][0]['source_volume_sense']='Both'
        with self.assertRaisesRegex(ValueError,'sense'):bundle.check_capture(mesh,self.native)

    def test_stale_receipt_stops_before_output_creation(self):
        # Test-only synthetic receipt, not a claim of executed native acquisition.
        fake=Path(self.tmp.name)/'synthetic-not-native.dll';fake.write_bytes(b'SYNTHETIC TEST ONLY; NOT A BINARY')
        mesh=copy.deepcopy(self.mesh);mesh.update(exporter_module_path=str(fake),cubit_version='SYNTHETIC-NO-EXECUTION')
        capture=Path(self.tmp.name)/'native-capture.json';capture.write_bytes(bundle.encoded(mesh))
        alternate=Path(self.tmp.name)/'alternate-mesh.e';alternate.write_bytes(self.native.read_bytes())
        receipt=Path(self.tmp.name)/'receipt.json';receipt.write_bytes(bundle.encoded({
            'complete':True,'batch_nographics':True,'capture_sha256':'0'*64,
            'alternate_sha256':bundle.digest(alternate.read_bytes()),'loaded_module_path':str(fake),
            'candidate_sha256':bundle.digest(fake.read_bytes()),'candidate_after_sha256':bundle.digest(fake.read_bytes()),
            'collector_source_sha256':bundle.digest(Path(__file__).read_bytes()),
            'native_source_digest':bundle.digest(b'synthetic source'),
            'source_journal_sha256':bundle.digest(b'synthetic journal'),
            'execution_id':'SYNTHETIC-CONTRACT-NO-EXECUTION','started_utc':'2000-01-01T00:00:00Z',
            'finished_utc':'2000-01-01T00:00:01Z'}))
        out=Path(self.tmp.name)/'output'
        with self.assertRaisesRegex(ValueError,'changed after receipt'):
            bundle.seal_capture(capture,receipt,out,schema_path=SOURCE.parents[2]/'schemas/mesh-bundle.schema.json')
        self.assertFalse(out.exists())

    def test_strict_json_rejects_duplicate_and_nonfinite(self):
        for value in (b'{"id":1,"id":2}',b'{"x":NaN}'):
            with self.assertRaises(ValueError):bundle.strict_json(value)

    def test_label_catalog_independent_exodus_membership_and_name(self):
        import numpy as np
        from scipy.io import netcdf_file
        mesh=copy.deepcopy(self.mesh)
        mesh['sets']={'sidesets':[{'id':4,'name':'tagged_face','kind':'face','edges':[],
                                  'faces':[mesh['faces'][0]],'node_ids':[1,2,3,4]}],
                      'nodesets':[{'id':8,'name':'tagged_nodes','node_ids':[1,2]}]}
        with netcdf_file(self.native,'a') as db:
            for key,size in [('len_name',32),('num_side_sets',1),('num_side_ss1',1),
                             ('num_node_sets',1),('num_nod_ns1',2)]:db.createDimension(key,size)
            def label(var,dim,text):
                db.createVariable(var,'c',(dim,'len_name'))[:] = np.array([list(text.ljust(32,'\0'))],dtype='S1')
            label('eb_names','num_el_blk','box_region')
            label('ss_names','num_side_sets','tagged_face')
            label('ns_names','num_node_sets','tagged_nodes')
            db.createVariable('ss_prop1','i',('num_side_sets',))[:] = [4]
            db.createVariable('elem_ss1','i',('num_side_ss1',))[:] = [1]
            db.createVariable('side_ss1','i',('num_side_ss1',))[:] = [5]
            db.createVariable('ns_prop1','i',('num_node_sets',))[:] = [8]
            db.createVariable('node_ns1','i',('num_nod_ns1',))[:] = [1,2]
        capture=Path(self.tmp.name)/'native-capture.json'
        alternate=Path(self.tmp.name)/'alternate-mesh.e';alternate.write_bytes(self.native.read_bytes())
        receipt=Path(self.tmp.name)/'receipt.json'
        def reseal():
            capture.write_bytes(bundle.encoded(mesh))
            receipt.write_bytes(bundle.encoded({'complete':True,'batch_nographics':True,
                'capture_sha256':bundle.digest(capture.read_bytes()),
                'alternate_sha256':bundle.digest(alternate.read_bytes())}))
        reseal()
        catalog=bundle.label_catalog(capture,receipt)
        self.assertEqual(catalog['sets'],mesh['sets'])
        self.assertTrue(catalog['native_label_membership_verified'])
        self.assertFalse(catalog['physical_conditions_transferred'])
        mesh['sets']['sidesets'][0]['name']='changed_name';reseal()
        with self.assertRaisesRegex(ValueError,'sideset labels'):bundle.label_catalog(capture,receipt)
        with self.assertRaisesRegex(ValueError,'cannot silently drop'):
            bundle.seal_capture(capture,receipt,Path(self.tmp.name)/'out',schema_path='not-read')


if __name__ == '__main__':
    unittest.main()
