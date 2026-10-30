"""The pre-analysis record must exist and match every release input before inference."""
from pathlib import Path
import json,hashlib

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        while b:=f.read(8388608):h.update(b)
    return h.hexdigest()

def require_freeze(root):
    root=Path(root);p=root/'manifests/freeze_record.json'
    if not p.exists():raise RuntimeError('Main inference blocked: reviewed prospective freeze record is absent.')
    rec=json.loads(p.read_text())
    private=rec.get('method')=='private_prospective_snapshot'
    if private:
        if not rec.get('user_authorized_private_completion') or not rec.get('approved_prospective_amendment') or rec.get('independent_external_timestamp') is not False:raise RuntimeError('Private freeze authorization/disclosure is incomplete')
    elif not rec.get('verified_external_timestamp') or not rec.get('user_authorized_destination'):raise RuntimeError('Freeze is unverified or unauthorized')
    if rec.get('protocol_version')!='0.4-draft':raise RuntimeError('Freeze version mismatch')
    # Any alternative timestamp method requires a prospective recorded amendment.
    if rec.get('method')!='public_exact_files' and not rec.get('approved_prospective_amendment'):raise RuntimeError('Alternative freeze rule has no approved amendment')
    if rec.get('method')=='public_exact_files' and not rec.get('unauthenticated_retrieval_verified'):raise RuntimeError('Public access unverified')
    m=root/rec['manifest_path']
    if sha256(m)!=rec['manifest_sha256']:raise RuntimeError('Frozen manifest changed')
    for item in json.loads(m.read_text())['files']:
        if sha256(root/item['path'])!=item['sha256']:raise RuntimeError('Frozen file changed: '+item['path'])
    return rec
