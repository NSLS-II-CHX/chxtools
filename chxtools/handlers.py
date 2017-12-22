# handler registration and database instantiation should be done
# here and only here!
from databroker import Broker
from databroker.assets.handlers_base import HandlerBase
#from chxtools.pims_readers.eiger import EigerImages
from eiger_io.fs_handler_dask import EigerHandlerDask, EigerImagesDask
from eiger_io.fs_handler import EigerHandler as EigerHandlerPIMS, EigerImages as EigerImagesPIMS

db = Broker.named('chx')

'''
Tried to allow function to change namespace did not work. 
DO NOT USE
'''
# toggle use of dask or no dask
# TODO : eventually choose one of the two
def use_pims():
    global EigerImages, EigerHandler, db
    EigerImages = EigerImagesPIMS
    EigerHandler = EigerHandlerPIMS
    db.reg.register_handler('AD_EIGER2', EigerHandler, overwrite=True)
    db.reg.register_handler('AD_EIGER', EigerHandler, overwrite=True)


def use_dask():
    global EigerImages, EigerHandler, db

    EigerImages = EigerImagesDask
    EigerHandler = EigerHandlerDask
    db.reg.register_handler('AD_EIGER2', EigerHandler, overwrite=True)
    db.reg.register_handler('AD_EIGER', EigerHandler, overwrite=True)

# call use_pims or use_dask
# default is use_dask()
# TODO : This is hard coded
# calling this after import won't change things, need to find a better way
use_pims()
