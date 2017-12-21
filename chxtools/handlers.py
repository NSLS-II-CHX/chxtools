# this does handler registration?
from databroker import Broker
from databroker.assets.handlers_base import HandlerBase
#from chxtools.pims_readers.eiger import EigerImages
from eiger_io.fs_handler_dask import EigerDaskHandler, EigerImagesDask
from eiger_io.fs_handler import EigerHandler, EigerImages

# toggle use of dask or no dask
use_dask = True
if use_dask:
    EigerHandler = EigerDaskHandler
    EigerImages = EigerImagesDask


db = Broker.named('chx')
db.reg.register_handler('AD_EIGER2', EigerHandler)
db.reg.register_handler('AD_EIGER2', EigerHandler)
