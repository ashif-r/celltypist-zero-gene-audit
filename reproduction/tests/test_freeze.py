import pytest,json
from src.freeze import require_freeze,sha256

def test_no_freeze_blocks_before_output(tmp_path):
 with pytest.raises(RuntimeError,match='absent'):require_freeze(tmp_path)
 assert not (tmp_path/'results').exists()
def test_changed_frozen_file_blocks(tmp_path):
 (tmp_path/'manifests').mkdir();p=tmp_path/'protocol.md';p.write_text('original')
 m=tmp_path/'manifests/files.json';m.write_text(json.dumps({'files':[{'path':'protocol.md','sha256':sha256(p)}]}))
 rec={'verified_external_timestamp':True,'user_authorized_destination':True,'protocol_version':'0.4-draft','method':'public_exact_files','unauthenticated_retrieval_verified':True,'manifest_path':'manifests/files.json','manifest_sha256':sha256(m)}
 (tmp_path/'manifests/freeze_record.json').write_text(json.dumps(rec));assert require_freeze(tmp_path)==rec
 p.write_text('changed')
 with pytest.raises(RuntimeError,match='changed'):require_freeze(tmp_path)
def test_private_commit_is_not_public_freeze(tmp_path):
 (tmp_path/'manifests').mkdir();(tmp_path/'manifests/freeze_record.json').write_text(json.dumps({'verified_external_timestamp':True,'user_authorized_destination':True,'protocol_version':'0.4-draft','method':'private_git_commit'}))
 with pytest.raises(RuntimeError,match='amendment'):require_freeze(tmp_path)


def test_private_prospective_amendment_is_honestly_disclosed(tmp_path):
 (tmp_path/'manifests').mkdir();p=tmp_path/'protocol.md';p.write_text('private rule')
 m=tmp_path/'manifests/files.json';m.write_text(json.dumps({'files':[{'path':'protocol.md','sha256':sha256(p)}]}))
 rec={'protocol_version':'0.4-draft','method':'private_prospective_snapshot','user_authorized_private_completion':True,'approved_prospective_amendment':True,'independent_external_timestamp':False,'manifest_path':'manifests/files.json','manifest_sha256':sha256(m)}
 (tmp_path/'manifests/freeze_record.json').write_text(json.dumps(rec));assert require_freeze(tmp_path)==rec
 rec['independent_external_timestamp']=True;(tmp_path/'manifests/freeze_record.json').write_text(json.dumps(rec))
 with pytest.raises(RuntimeError,match='disclosure'):require_freeze(tmp_path)
