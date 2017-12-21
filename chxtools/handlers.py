# this does handler registration?
from databroker import Broker
from databroker.assets.handlers_base import HandlerBase
#from chxtools.pims_readers.eiger import EigerImages
from eiger_io.fs_handlers_dask import EigerDaskHandler
from eiger_io.fs_handlers import EigerHandler

# toggle use of dask or no dask
EigerHandler = EigerDaskHandler


db = Broker.named('chx')
db.reg.register_handler('AD_EIGER2', EigerHandler)
db.reg.register_handler('AD_EIGER2', EigerHandler)
