# handler registration and database instantiation should be done
# here and only here!
from databroker import Broker
from databroker.assets.handlers_base import HandlerBase
#from chxtools.pims_readers.eiger import EigerImages
from eiger_io.fs_handler_dask import EigerHandlerDask, EigerImagesDask
from eiger_io.fs_handler import EigerHandler, EigerImages

# toggle use of dask or no dask
# TODO : eventually choose one of the two
use_dask = False
if use_dask:
    EigerHandler = EigerHandlerDask
    EigerImages = EigerImagesDask


db = Broker.named('chx')
db.reg.register_handler('AD_EIGER2', EigerHandler)
db.reg.register_handler('AD_EIGER2', EigerHandler)
