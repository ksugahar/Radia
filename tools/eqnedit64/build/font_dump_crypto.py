"""Encrypt private CI crash evidence before uploading to a public repository."""
import argparse
import base64
import hashlib
import io
import os
from pathlib import Path
import struct
import zipfile

from cryptography.hazmat.primitives import hashes, serialization
from cryptography.hazmat.primitives.asymmetric import padding, rsa
from cryptography.hazmat.primitives.ciphers.aead import AESGCM

MAGIC = b'EQN-DUMP-1\0'


def oaep():
    return padding.OAEP(mgf=padding.MGF1(hashes.SHA256()), algorithm=hashes.SHA256(), label=MAGIC)


def encrypt(public, data):
    key = AESGCM.generate_key(bit_length=256)
    nonce = os.urandom(12)
    wrapped = public.encrypt(key, oaep())
    return MAGIC + struct.pack('>H', len(wrapped)) + wrapped + nonce + AESGCM(key).encrypt(nonce, data, MAGIC)


def decrypt(private, blob):
    if not blob.startswith(MAGIC):
        raise ValueError('Unknown evidence envelope')
    start = len(MAGIC) + 2
    size = struct.unpack('>H', blob[len(MAGIC):start])[0]
    key = private.decrypt(blob[start:start+size], oaep())
    nonce = blob[start+size:start+size+12]
    return AESGCM(key).decrypt(nonce, blob[start+size+12:], MAGIC)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('command', choices=['generate', 'encrypt', 'decrypt', 'selftest'])
    parser.add_argument('--path', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--key', type=Path)
    args = parser.parse_args()
    if args.command in ('generate', 'selftest'):
        private = rsa.generate_private_key(public_exponent=65537, key_size=3072)
        if args.command == 'selftest':
            data = b'private diagnostic fixture'
            blob = encrypt(private.public_key(), data)
            assert decrypt(private, blob) == data
            try:
                decrypt(private, blob[:-1] + bytes([blob[-1] ^ 1]))
            except Exception:
                pass
            else:
                raise AssertionError('Tampering accepted')
            print('Round-trip and tamper rejection PASS')
            return
        with args.path.open('xb') as stream:
            stream.write(private.private_bytes(serialization.Encoding.PEM, serialization.PrivateFormat.PKCS8, serialization.NoEncryption()))
        public = private.public_key().public_bytes(serialization.Encoding.DER, serialization.PublicFormat.SubjectPublicKeyInfo)
        args.path.with_suffix('.public.txt').write_text(base64.b64encode(public).decode('ascii'), encoding='ascii')
        print('Created private key locally; public fingerprint:', hashlib.sha256(public).hexdigest())
    elif args.command == 'encrypt':
        public = serialization.load_der_public_key(base64.b64decode(os.environ['FONT_DUMP_PUBLIC_KEY'], validate=True))
        files = sorted(p for p in args.path.iterdir() if p.is_file())
        if not files or sum(p.stat().st_size for p in files) > 1024**3:
            raise ValueError('Missing or oversized evidence')
        archive = io.BytesIO()
        with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as zipped:
            for path in files:
                zipped.write(path, path.name)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open('xb') as stream:
            stream.write(encrypt(public, archive.getvalue()))
        print('Encrypted', len(files), 'files; plaintext is not uploaded.')
    else:
        private = serialization.load_pem_private_key(args.key.read_bytes(), password=None)
        with args.output.open('xb') as stream:
            stream.write(decrypt(private, args.path.read_bytes()))
        print('Decrypted locally to', args.output)


if __name__ == '__main__':
    main()
