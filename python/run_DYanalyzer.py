import FWCore.ParameterSet.Config as cms

process = cms.Process("DYLHE")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )



#process.TFileService = cms.Service("TFileService",
                                   #fileName = cms.string("/home/jhchoi/mcm_filter_test/HIG-RunIIFall18wmLHEGS-01691.root"),
#                                   closeFileFast = cms.untracked.bool(True)
#                                   )
        
process.source = cms.Source("PoolSource",
                                    # replace 'myfile.root' with the source file you want to use                                                                      
           fileNames = cms.untracked.vstring(

"file:/home/jhchoi/mcm_filter_test/HIG-RunIIFall18wmLHEGS-01691_inLHE.root"
      )
)
       
process.DYana = cms.EDAnalyzer('DYanalyzer',
                              
                              genSrc = cms.InputTag("genParticles")
#prunedGenParticles
)


process.p = cms.Path(process.DYana)
