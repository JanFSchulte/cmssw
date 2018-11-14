import FWCore.ParameterSet.Config as cms

from RecoTracker.TkNavigation.NavigationSchoolESProducer_cfi import *
telescopeNavigationSchoolESProducer = navigationSchoolESProducer.clone()
telescopeNavigationSchoolESProducer.ComponentName       = cms.string('TelescopeNavigationSchool')
telescopeNavigationSchoolESProducer.SimpleMagneticField = cms.string('')


