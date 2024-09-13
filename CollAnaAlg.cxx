#include "CollAnaAlg/CollAnaAlg.h"
#include "CollAnaAlg/TrackMatch.h"
#include "McDecayModeSvc/McDecayModeSvc.h"
#include "BestDTagSvc/BestDTagSvc.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "EvtRecEvent/EvtRecDTag.h"
#include "EvtRecEvent/EvtRecVeeVertex.h"
#include "EvtRecEvent/EvtRecPi0.h"
#include "DstEvent/TofHitStatus.h"

#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Geometry/Point3D.h"


#include "McTruth/McParticle.h"
#include "EvTimeEvent/RecEsTime.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;

//#include "DTagTool/DTagTool.h"
#include "VertexFit/Helix.h" 
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "ParticleID/ParticleID.h"
#include "CLHEP/Geometry/Point3D.h"
#include "SimplePIDSvc/ISimplePIDSvc.h"
#include "McDecayModeSvc/McDecayModeSvc.h"
#include "GaudiKernel/IPartPropSvc.h"

#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

#include <vector>
const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};
typedef std::vector<int> Vint;
typedef std::vector<double> Vdouble;
typedef std::vector<HepLorentzVector> Vp4;

static long m_cout_all(0), m_cout_dcol(0), m_cout_pass(0);


CollAnaAlg::CollAnaAlg(const std::string& name, ISvcLocator* pSvcLocator):Algorithm(name,pSvcLocator){
	declareProperty("Vr0cut",   m_vr0cut=1.0);
	declareProperty("Vz0cut",   m_vz0cut=10.0);
	declareProperty("Dang_cut",m_dang_cut=20);
	declareProperty("TrkCos_cut",m_trk_cos_cut=0.93);
	declareProperty("Energy_cut_b",m_energy_cut_b=0.025);
	declareProperty("Energy_cut_e",m_energy_cut_e=0.05);
	declareProperty("Debug",m_debug=0);
	declareProperty("Nele_cut",m_nele_cut=0); /// charged pion
	declareProperty("Num_cut",numParticle_cut=3);
	declareProperty("ApplyBoost", m_applyBoost = true);
	declareProperty("ApplyTimeCut",m_applyTimeCut=true);
	declareProperty("ApplyDangCut",m_applyDangCut=true);
	declareProperty("beamE",m_beamE=2.13);
	declareProperty("ReadBeamEFromDB",m_ReadBeamEFromDB=true);
	declareProperty("UseCalibBeamE",m_usecalibBeamE=false);
	declareProperty("UseTruthAxis",m_UseTruthAxis=false);
	declareProperty("Isqqbar",     m_isqqbar = true);

}

CollAnaAlg::~CollAnaAlg(){
	//add your code for deconstructor
}

StatusCode CollAnaAlg::initialize(){
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"CollAnaAlg::initialize()"<<endreq;

	m_irun=-100;

	//add your code here
	StatusCode status;
	NTuplePtr nt2(ntupleSvc(), "FILE1/mctruth");
	if ( nt2 ) m_tuple2 = nt2;
	else {
		m_tuple2 = ntupleSvc()->book ("FILE1/mctruth", CLID_ColumnWiseTuple, "exam N-Tuple example");
		if ( m_tuple2 )    {
			status = m_tuple2->addItem ("run",   m_run_2);
			status = m_tuple2->addItem ("event", m_event_2);
			status = m_tuple2->addItem ("Nleaf", m_nleaf_2, 0, 35);
			status = m_tuple2->addItem ("gluon", m_gluon_2);
			status = m_tuple2->addItem ("q_p4", 4, m_p4_q_2);
			status = m_tuple2->addItem ("qbar_p4",4, m_p4_qbar_2);
			status = m_tuple2->addItem ("qid", m_qid_2);
			status = m_tuple2->addItem ("qbarid", m_qbarid_2);
			status = m_tuple2->addItem("ntotMC",   m_tot_2, 0,100);
			status = m_tuple2->addIndexedItem("pid",   m_tot_2, m_pid_2);
			status = m_tuple2->addIndexedItem("mother",m_tot_2, m_mother_2);
			status = m_tuple2->addIndexedItem("p4mc",m_tot_2, 4, m_p4mc_2);
			status = m_tuple2->addItem("Nleaf_charged",     m_nleaf_ch_2, 0, 20);
			status = m_tuple2->addItem("Nleaf_ne",     m_nleaf_ne_2);
			status = m_tuple2->addItem("Npi0_mc",     m_n_pi0_2);
			status = m_tuple2->addItem("npion_mc",     m_npion_mc_2);
			status = m_tuple2->addItem("net_leaf_charge",     m_net_leaf_charge_2);

			status = m_tuple2->addIndexedItem("Leaf_info",  m_nleaf_ch_2, 7, m_mom_leaf_2);//fourmom of leaf particles and its charge and pid and its mother
			status = m_tuple2->addIndexedItem("Leaf_cosTheta",  m_nleaf_ch_2,  m_mc_cos_2);// costheta of a charged leaf particle //befor boost
			status = m_tuple2->addIndexedItem("Leaf_phi",  m_nleaf_ch_2,  m_mc_phi_2);// theta of a charged leaf particle //befor boost
			status = m_tuple2->addIndexedItem("Leaf_theta",  m_nleaf_ch_2,  m_mc_theta_2);// costheta of a charged leaf particle
			status = m_tuple2->addIndexedItem ("belong_mc", m_nleaf_ch_2, m_belong_mc_2);//belong to q = 1 , belong to qbar=-1
			status = m_tuple2->addIndexedItem ("charge_mc", m_nleaf_ch_2, m_charge_mc_2);//charge of the mc particle
			status = m_tuple2->addIndexedItem ("z_mc",      m_nleaf_ch_2, m_z_mc_2);//z of the mc

			status = m_tuple2->addIndexedItem ("angle_mc",  m_nleaf_ch_2,20 ,m_angle_mc_2);//angle between any given pi-pair // after boost
			status = m_tuple2->addIndexedItem ("isfavor_mc", m_nleaf_ch_2,20 ,m_isfavor_mc_2);
			status = m_tuple2->addIndexedItem ("phi0_mc",  m_nleaf_ch_2,20 ,m_phi0_mc_2);//2*phi0 between any given pi-pair
			status = m_tuple2->addIndexedItem ("mass_mc",  m_nleaf_ch_2,20 ,m_mass_mc_2);//possible resonsance
			status = m_tuple2->addIndexedItem ("qt_mc",     m_nleaf_ch_2,20, m_qt_mc_2);//qt of the p1 vs p2 
		}
		else    {
			log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple2) << endmsg;
			return StatusCode::FAILURE;
		}
	}

	NTuplePtr nt1(ntupleSvc(), "FILE1/qq");
	if ( nt1 ) m_tuple1 = nt1;
	else {
		m_tuple1 = ntupleSvc()->book ("FILE1/qq", CLID_ColumnWiseTuple, "exam N-Tuple example");
		if ( m_tuple1 )    {
			status = m_tuple1->addItem ("run",   m_run);
			status = m_tuple1->addItem ("event", m_event);
			status = m_tuple1->addItem ("Ebeam", m_Ebeam);
			status = m_tuple1->addItem ("index", m_index);

			status = m_tuple1->addItem("Nleaf",          m_nleaf, 0, 40);
			status = m_tuple1->addIndexedItem("leafId",  m_nleaf, m_leafId);
			status = m_tuple1->addIndexedItem("p4_leaf",  m_nleaf, 4, m_p4_leaf);
			//infomation from charged lead pariticles
			status = m_tuple1->addItem("Nleaf_charged",     m_nleaf_ch, 0, 20);
			status = m_tuple1->addItem("Nleaf_ne",     m_nleaf_ne); //netural particles;
			status = m_tuple1->addItem("Npi0_mc",     m_n_pi0);
			status = m_tuple1->addItem("npion_mc",     m_npion_mc);
			//status = m_tuple1->addItem("net_leaf_charge",     m_net_leaf_charge);
			status = m_tuple1->addIndexedItem("Leaf_cosTheta",  m_nleaf_ch,  m_mc_cos);// costheta of a charged leaf particle befor boost(in lab frame)
			status = m_tuple1->addIndexedItem("Leaf_theta",  m_nleaf_ch,  m_mc_theta);// theta of a charged leaf particle befor boost(in lab frame)
			status = m_tuple1->addIndexedItem("Leaf_phi",   m_nleaf_ch,  m_mc_phi);// phi of a charged leaf particle befor boost(in lab frame)
			status = m_tuple1->addIndexedItem("Leaf_info",  m_nleaf_ch, 7, m_mom_leaf);//fourmom of leaf particles (in cms sys) and its charge and pid and its mother
			status = m_tuple1->addIndexedItem ("belong_mc", m_nleaf_ch, m_belong_mc);//belong to q = 1 , belong to qbar=-1
			status = m_tuple1->addIndexedItem ("trkIndex_mc",m_nleaf_ch,m_trkIndex_mc);//used to identify the track
			status = m_tuple1->addIndexedItem ("charge_mc", m_nleaf_ch, m_charge_mc);//charge of the mc particle
			status = m_tuple1->addIndexedItem ("z_mc",      m_nleaf_ch, m_z_mc);//z of the mc
			status = m_tuple1->addIndexedItem ("pt_mc",     m_nleaf_ch, m_pt_mc);//pt to the Truth thrust axis
			status = m_tuple1->addIndexedItem ("phi_mc",    m_nleaf_ch, m_phi_mc);//angle between Thrust and any pi 
			status = m_tuple1->addIndexedItem ("Ifrec",     m_nleaf_ch, m_Ifrec); 

			status = m_tuple1->addIndexedItem ("Qt_mc",     m_nleaf_ch,20, m_qt_mc);//pt to the Truth thrust axis
			status = m_tuple1->addIndexedItem ("angle_mc",  m_nleaf_ch,20 ,m_angle_mc);//angle between any given pi-pair // // after boost
			status = m_tuple1->addIndexedItem ("isfavor_mc", m_nleaf_ch,20 ,m_isfavor_mc);
			status = m_tuple1->addIndexedItem ("phi0_mc",  m_nleaf_ch,20 ,m_phi0_mc);//2*phi0 between any given pi-pair
			status = m_tuple1->addIndexedItem ("mass_mc",  m_nleaf_ch,20 ,m_mass_mc);//2*phi0 between any given pi-pair

			status = m_tuple1->addItem("nMC",              m_idxmc, 0, 100);
			status = m_tuple1->addIndexedItem("pdgid",     m_idxmc, m_pdgid);
			status = m_tuple1->addIndexedItem("motheridx", m_idxmc, m_motheridx);

			status = m_tuple1->addItem("p4_q",    4, m_p4_q);
			status = m_tuple1->addItem("p4_qbar", 4, m_p4_qbar);
			//status = m_tuple1->addItem("thru_mc",  m_thru_mc);
			status = m_tuple1->addItem ("gluon", m_gluon);//number of gluon
			status = m_tuple1->addItem ("qid", m_qid);
			status = m_tuple1->addItem ("qbarid", m_qbarid);
			status = m_tuple1->addItem("ntotMC",   m_tot, 0, 100);
			status = m_tuple1->addIndexedItem("pid",   m_tot, m_pid);
			status = m_tuple1->addIndexedItem("mother",m_tot, m_mother);
			status = m_tuple1->addIndexedItem("motherindex", m_tot, m_motherindex);
			status = m_tuple1->addIndexedItem("p4mc", m_tot, 4, m_p4mc);
			status = m_tuple1->addIndexedItem("mc_primaryParticle", m_tot, m_mc_primaryParticle);
			status = m_tuple1->addIndexedItem("mc_leafParticle", m_tot, m_mc_leafParticle);
			status = m_tuple1->addIndexedItem("mc_decayFromGen", m_tot, m_mc_decayFromGen);
			status = m_tuple1->addIndexedItem("mc_decayInFlight", m_tot, m_mc_decayInFlight);
			status = m_tuple1->addIndexedItem("mc_children", m_tot, m_mc_children);
			status = m_tuple1->addIndexedItem("mc_childrenFlight", m_tot, m_mc_childrenFlight);
			status = m_tuple1->addIndexedItem("mc_isFinal", m_tot, m_mc_isFinal);

			status = m_tuple1->addItem("Thr_diy_mc", m_Thr_diy_mc);
			status = m_tuple1->addItem("Thr_axis_mc",3,m_Thr_axis_mc);

			status = m_tuple1->addItem("ngood", m_ngood,0,20);
			status = m_tuple1->addItem("npion", m_npion,0,20);
			status = m_tuple1->addItem("nkaon", m_nkaon,0,10);
			status = m_tuple1->addItem("nele", m_nele);
			status = m_tuple1->addItem("nproton", m_nproton);
			status = m_tuple1->addItem("nshower", m_nshower,0,30);

			status = m_tuple1->addIndexedItem ("P4_ch", m_ngood, 4, m_fourmom_ch);//four momentum after boost
			status = m_tuple1->addIndexedItem ("P_ch",  m_ngood, m_mom_ch);//after boost

			status = m_tuple1->addIndexedItem ("pion_cosT",       m_npion,     m_pion_cosT); //befor boost
			status = m_tuple1->addIndexedItem ("pion_Theta",      m_npion,     m_pion_Theta);//befor boost
			status = m_tuple1->addIndexedItem ("pion_Phi",        m_npion,     m_pion_phi);// befor boost
			status = m_tuple1->addIndexedItem ("pion_belong",     m_npion,     m_pion_belong); //belong Thrust +or-
			status = m_tuple1->addIndexedItem ("pion_charge",     m_npion,     m_pion_charge);
			status = m_tuple1->addIndexedItem ("pion_theta0",     m_npion,     m_pion_theta0); //theta after boost
			status = m_tuple1->addIndexedItem ("pion_pt",         m_npion,     m_pion_pt);//pt after boost
			status = m_tuple1->addIndexedItem ("pion_trkIndex",   m_npion,     m_pion_trkIndex);//trkIndex from match
			status = m_tuple1->addIndexedItem ("pion_type_match", m_npion,     m_pion_type_match);//particle type from truth match
			status = m_tuple1->addIndexedItem ("pion_z",          m_npion,     m_pion_z);
			status = m_tuple1->addIndexedItem ("pion_phi1",       m_npion ,    m_pion_phi1);// phi of each pion
			status = m_tuple1->addIndexedItem ("pion_p4",         m_npion, 4,  m_pion_p4);
			status = m_tuple1->addIndexedItem ("pion_angle",      m_npion,20,  m_pion_angle);
			status = m_tuple1->addIndexedItem ("pion_mass",       m_npion,20,  m_pion_mass);//// possible resonsance 
			status = m_tuple1->addIndexedItem ("pion_isfavor",    m_npion,20,  m_pion_isfavor);
			status = m_tuple1->addIndexedItem ("pion_ptvs",       m_npion,20,  m_pion_ptvs);
			status = m_tuple1->addIndexedItem ("pion_Qt",         m_npion,20,  m_pion_qtvs);
			status = m_tuple1->addIndexedItem ("pion_phi0",       m_npion,20,  m_pion_phi0);//// 2phi0 of each pion pair

			status = m_tuple1->addIndexedItem ("kaon_cosT",       m_nkaon,     m_kaon_cosT); //befor boost
			status = m_tuple1->addIndexedItem ("kaon_Theta",      m_nkaon,     m_kaon_Theta);//befor boost
			status = m_tuple1->addIndexedItem ("kaon_Phi",        m_nkaon,     m_kaon_phi);// befor boost
			status = m_tuple1->addIndexedItem ("kaon_belong",     m_nkaon,     m_kaon_belong); //belong Thrust +or-
			status = m_tuple1->addIndexedItem ("kaon_charge",     m_nkaon,     m_kaon_charge);
			status = m_tuple1->addIndexedItem ("kaon_theta0",     m_nkaon,     m_kaon_theta0); //theta after boost
			status = m_tuple1->addIndexedItem ("kaon_pt",         m_nkaon,     m_kaon_pt);//pt after boost
			status = m_tuple1->addIndexedItem ("kaon_trkIndex",   m_nkaon,     m_kaon_trkIndex);//trkIndex from match
			status = m_tuple1->addIndexedItem ("kaon_type_match", m_nkaon,     m_kaon_type_match);//particle type from truth match
			status = m_tuple1->addIndexedItem ("kaon_z",          m_nkaon,     m_kaon_z);
			status = m_tuple1->addIndexedItem ("kaon_phi1",       m_nkaon ,    m_kaon_phi1);// phi of each kaon
			status = m_tuple1->addIndexedItem ("kaon_p4",         m_nkaon, 4,  m_kaon_p4);
			status = m_tuple1->addIndexedItem ("kaon_angle",      m_nkaon,10,  m_kaon_angle);
			status = m_tuple1->addIndexedItem ("kaon_mass",       m_nkaon,10,  m_kaon_mass);//// possible resonsance 
			status = m_tuple1->addIndexedItem ("kaon_isfavor",    m_nkaon,10,  m_kaon_isfavor);
			status = m_tuple1->addIndexedItem ("kaon_ptvs",       m_nkaon,10,  m_kaon_ptvs);
			status = m_tuple1->addIndexedItem ("kaon_Qt",         m_nkaon,10,  m_kaon_qtvs);
			status = m_tuple1->addIndexedItem ("kaon_phi0",       m_nkaon,10,  m_kaon_phi0);//// 2phi0 of each kaon pair

			status = m_tuple1->addIndexedItem ("kpi_angle",       m_nkaon,20,  m_kpi_angle);
			status = m_tuple1->addIndexedItem ("kpi_mass",        m_nkaon,20,  m_kpi_mass);//// possible resonsance 
			status = m_tuple1->addIndexedItem ("kpi_isfavor",     m_nkaon,20,  m_kpi_isfavor);
			status = m_tuple1->addIndexedItem ("kpi_ptvs",        m_nkaon,20,  m_kpi_ptvs);
			status = m_tuple1->addIndexedItem ("kpi_Qt",          m_nkaon,20,  m_kpi_qtvs);
			status = m_tuple1->addIndexedItem ("kpi_phi0",        m_nkaon,20,  m_kpi_phi0);//// 2phi0 of each kpi pair
			status = m_tuple1->addIndexedItem ("pik_ptvs",        m_npion,10,  m_pik_ptvs);
			status = m_tuple1->addIndexedItem ("pik_Qt",          m_npion,10,  m_pik_qtvs);
			status = m_tuple1->addIndexedItem ("pik_phi0",        m_npion,10,  m_pik_phi0);//// 2phi0 of each pik pair


			status = m_tuple1->addIndexedItem ("P4_ne",  m_nshower, 4, m_fourmom_ne);
			status = m_tuple1->addItem("Evis",  m_evis);

			//separate pion1 and pion2
			status = m_tuple1->addItem("Npion1", m_npion1);
			status = m_tuple1->addItem("Npion2", m_npion2);
			status = m_tuple1->addItem("Nkaon1", m_nkaon1);
			status = m_tuple1->addItem("Nkaon2", m_nkaon2);
			//status = m_tuple1->addItem("Npion1", m_npion1,0,20);
			//status = m_tuple1->addItem("Npion2", m_npion2,0,20);
			//status = m_tuple1->addItem("Nkaon1", m_nkaon1,0,20);
			//status = m_tuple1->addItem("Nkaon2", m_nkaon2,0,20);

			status = m_tuple1->addItem("Thr_diy", m_Thr_diy);
			status = m_tuple1->addItem("Thr_axis",3,m_Thr_axis);
			status = m_tuple1->addItem("theta",m_Thr_theta);///theta of Thr_axis

		}
		else    {
			log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple1) << endmsg;
			return StatusCode::FAILURE;
		}
	}

	NTuplePtr nt3(ntupleSvc(), "FILE1/check");
	if ( nt3 ) m_tuple3 = nt3;
	else {
		m_tuple3 = ntupleSvc()->book ("FILE1/check", CLID_ColumnWiseTuple, "exam N-Tuple example");
		if ( m_tuple3 )    {
			status = m_tuple3->addItem ("run",   m_run_3);
			status = m_tuple3->addItem ("event", m_event_3);
			status = m_tuple3->addItem ("gluon",  m_gluon_3, 0, 100 );
			status = m_tuple3->addItem ("nGent",  m_n_Gent4);
			status = m_tuple3->addIndexedItem("Egluon", m_gluon_3, m_Egluon);
			status = m_tuple3->addItem("TotP4", 4, m_totP4);
			status = m_tuple3->addItem("Z0", 4, m_Z0);
			status = m_tuple3->addItem("rISR", 4, m_rISR);
			status = m_tuple3->addItem("string", 4, m_string);
			status = m_tuple3->addItem("q_p4", 4, m_q_p4);
			status = m_tuple3->addItem("qbar_p4", 4, m_qbar_p4);
			status = m_tuple3->addItem("daughter", 4, m_daughter);
			status = m_tuple3->addItem("Fid", m_fid);
		}
		else    {
			log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple3) << endmsg;
			return StatusCode::FAILURE;
		}
	}

	if(m_UseTruthAxis) cout<<" info: UseTruthAxis is ture"<<endl;

	IPartPropSvc* p_PartPropSvc;
	static const bool CREATEIFNOTTHERE(true);
	StatusCode PartPropStatus = service("PartPropSvc", p_PartPropSvc, CREATEIFNOTTHERE);
	if (!PartPropStatus.isSuccess() || 0 == p_PartPropSvc) {
		log << MSG::ERROR << " Could not initialize Particle Properties Service" << endreq;
		return PartPropStatus;
	}

	m_particleTable = p_PartPropSvc->PDT();
	return StatusCode::SUCCESS;
}

StatusCode CollAnaAlg::beginRun(){
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"CollAnaAlg::beginRun()"<<endreq;
	//add your code here
	return StatusCode::SUCCESS;
}

StatusCode CollAnaAlg::execute(){
	m_cout_all++;
	m_index = m_cout_all;
	if(m_debug)cout<<"***************************Collins V32 begin execute event: "<<m_index<<"***************************"<<endl;
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"CollAnaAlg::execute()"<<endreq;

	Hep3Vector epem(0,0,1);
	Hep3Vector plane_mc(0,0,0);
	Hep3Vector plane_ref_mc(0,0,0);

	//initialize
	m_Thr_theta =-10;
	m_npion1 = 0;
	m_npion2 = 0;
	m_nkaon1 = 0;
	m_nkaon2 = 0;
	for(int ll =0 ;ll <4;ll++) {m_Z0[ll]=0;}
	for(int ll =0 ;ll <4;ll++) {m_totP4[ll]=0;  m_daughter[ll] =0 ; m_rISR[ll]=0;}
	for(int i = 0 ; i<20;i++) {
		for(int ll =0 ;ll <7;ll++) {m_mom_leaf[i][ll]=-10;m_mom_leaf_2[i][ll]=-10;}
		m_phi_mc[i] = -10;
		m_Ifrec[i] = -10;
		for(int j = 0 ; j<20;j++){
			m_phi0_mc[i][j] = -10;
			m_phi0_mc_2[i][j] = -10;
			m_mass_mc_2[i][j]= -10;
			m_mass_mc[i][j]= -10;
			m_qt_mc_2[i][j] = -10;
			m_qt_mc[i][j] = -10;
			m_isfavor_mc[i][j] =-10;
			m_isfavor_mc_2[i][j] =-10;
			m_angle_mc[i][j] = -10;
			m_angle_mc_2[i][j] = -10;
			m_mc_phi[i] = -10;
			m_mc_phi_2[i] = -10;
			m_mc_cos[i] = -10;
			m_mc_cos_2[i] = -10;
			m_mc_theta[i] = -10;
			m_mc_theta_2[i] = -10;
		}
	}
	for(int i = 0 ; i<20;i++) {
		m_pion_cosT[i] =-10; //data
		m_pion_Theta[i] =-10;//data
		m_pion_phi[i] =-10; //data
		m_pion_belong[i]=-10;
		m_pion_theta0[i]=-10;//data
		m_pion_pt[i]=-10;
		m_pion_trkIndex[i]=-10;
		m_pion_type_match[i]=-999;
		m_pion_z[i] =-10;//data
		m_pion_phi1[i] = -10; //data
		for (int j = 0; j<4; j++){
			m_pion_p4[i][j] = -10; //data
		}
		for(int j = 0 ; j<20;j++){
			m_pion_phi0[i][j] = -10;//data
			m_pion_mass[i][j] = -10;//data
			m_pion_isfavor[i][j] = -10; //data
			m_pion_angle[i][j]=-10; //data
			m_pion_ptvs[i][j] = -10; //data
			m_pion_qtvs[i][j] = -10; //data
		}
	}
	for(int i = 0 ; i<10;i++) {
		m_kaon_cosT[i] =-10; //data
		m_kaon_Theta[i] =-10;//data
		m_kaon_phi[i] =-10; //data
		m_kaon_belong[i]=-10;
		m_kaon_theta0[i]=-10;//data
		m_kaon_pt[i]=-10;
		m_kaon_trkIndex[i]=-10;
		m_kaon_type_match[i]=-999;
		m_kaon_z[i] =-10;//data
		m_kaon_phi1[i] = -10; //data
		for (int j = 0; j<4; j++){
			m_kaon_p4[i][j] = -10; //data
		}
		for(int j = 0 ; j<10;j++){
			m_kaon_phi0[i][j] = -10;//data
			m_kaon_mass[i][j] = -10;//data
			m_kaon_isfavor[i][j] = -10; //data
			m_kaon_angle[i][j]=-10; //data
			m_kaon_ptvs[i][j] = -10; //data
			m_kaon_qtvs[i][j] = -10; //data
		}
	}
	// kpi
	for (int i = 0; i< 10; i++){
		for (int j = 0; j< 20; j++){
			m_kpi_angle[i][j] = -10;
			m_kpi_mass[i][j] = -10;
			m_kpi_isfavor[i][j] = -10;
			m_kpi_ptvs[i][j] = -10;
			m_kpi_qtvs[i][j] = -10;
			m_kpi_phi0[i][j] = -10;
		}
	}
	//pik
	for (int i = 0; i< 20; i++){
		for (int j = 0; j< 10; j++){
			m_pik_ptvs[i][j] = -10;
			m_pik_qtvs[i][j] = -10;
			m_pik_phi0[i][j] = -10;
		}
	}

	if(m_debug) cout<<"read EventHeader"<<endl;
	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
	if(!eventHeader){
		log << MSG::ERROR << "EventHeader not found" << endreq;
		return  StatusCode::SUCCESS;
	}

	int runNo = eventHeader->runNumber();
	int eventNo = eventHeader->eventNumber();
	//if(m_debug) cout<<" runNo : "<< runNo <<", eventNo "<< eventNo<<endl;
	//cout<<" runNo : "<< runNo <<", eventNo "<< eventNo<<endl;
	//if(runNo!=9613||eventNo!=587787) return  StatusCode::SUCCESS;

	m_run = runNo;
	m_event = eventNo;
	m_run_2 = runNo;
	m_event_2 = eventNo;
	m_run_3 = runNo;
	m_event_3 = eventNo;

	if(m_debug) cout<<"read Ebeam from Database"<<endl;

	if(m_ReadBeamEFromDB && m_irun!=runNo){
		m_irun=runNo;
		if(m_usecalibBeamE) m_readDb.setcalib(true);
		m_beamE=m_readDb.getbeamE(m_irun,m_beamE);
	}
	m_Ebeam = m_beamE;
	if(m_debug) cout<<" m_Ebeam : "<< m_beamE<<endl;

	int nLeaf=0;
	int nLeaf_charged=0;
	int cnt_pi=0;
	int nLeaf_netural = 0;
	int nmc_pi0 = 0;

	Hep3Vector  Thr_axis_truth(0,0,0);

	if( m_cout_all % 10000 ==0 ) cout<<"event cont "<< m_cout_all<<endl;

	log << MSG::DEBUG <<"run, evtnum = "
		<< m_run << " , "
		<< m_event <<endreq;

	int MCmode=-99;
	std::vector<int> pdgid;
	std::vector<int> TrackIndex;
	std::vector<int> mother;
	std::vector<int> motherindex;
	pdgid.clear();
	TrackIndex.clear();
	mother.clear();
	motherindex.clear();

	IMcDecayModeSvc* i_svc;
	StatusCode sc_DecayModeSvc = service("McDecayModeSvc", i_svc);
	if ( sc_DecayModeSvc.isFailure() ){
		log << MSG::FATAL << "Could not load McDecayModeSvc!" << endreq;
		return sc_DecayModeSvc;
	}
	m_svc = dynamic_cast<McDecayModeSvc*>(i_svc);

	Hep3Vector  Thr_axis_mc(-1,-1,-1);
	const int Nbin_1=90;
	const int Nbin_2=180;
	double step_1= 0.5*CLHEP::pi/Nbin_1;
	double step_2= 2*(CLHEP::pi)/Nbin_2;
	double Thrust_diy_mc=0;
	double Thrust_record_mc[Nbin_1][Nbin_2];
	double Thrust_unit_mc[Nbin_1][Nbin_2];

	for(int aa=0;aa<Nbin_1;aa++){
		for(int bb=0;bb<Nbin_2;bb++){
			Thrust_record_mc[aa][bb] = 0; Thrust_unit_mc[aa][bb]=0;
		}
	}

	m_gluon=0;
	m_gluon_2=0;
	m_gluon_3=0;
	int ntot_MC=0;
	int n_Gent4=0;

	Vp4 mc_charged_hadron; mc_charged_hadron.clear();
	Vint mc_charge; mc_charge.clear();
	Vint mc_pid; mc_pid.clear();
	Vint mc_mother; mc_mother.clear();
	Vint mc_charged_trkIndex; mc_charged_trkIndex.clear();
	Vint mc_charged_belong; mc_charged_belong.clear();

	if (eventHeader->runNumber()<0){
		HepLorentzVector p4_q(0,0,0,0);
		HepLorentzVector p4_qbar(0,0,0,0);

		int cnt_q=0,cnt_qbar=0;
		int qid=0,qbarid=0;
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");

		int m_numParticle = 0;
		if (!mcParticleCol){
			std::cout << "Could not retrieve McParticelCol" << std::endl;
			return StatusCode::FAILURE;
		} else {
			Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
			for (; iter_mc != mcParticleCol->end(); iter_mc++) {
				if(m_debug){
					cout<<"ntot_MC "<<setw(2)<<ntot_MC<<" ";
					cout<<"MCParticleCol::";
					cout<<"isPrimaryPar("<<(*iter_mc)->primaryParticle()<<") ";
					cout<<"ParProperty("<<(*iter_mc)->particleProperty()<<") ";
					cout<<"decayFromGenerator("<<(*iter_mc)->decayFromGenerator()<<") ";
					cout<<"mother("<<(*iter_mc)->mother().particleProperty()<<") ";
					cout<<"motherindex("<<(*iter_mc)->mother().trackIndex()<<") ";
					cout<<"leafPar("<<(*iter_mc)->leafParticle()<<")"<<endl;
					//continue;
				}

				//ISR
				//Z , and ISR primaryParticle==1, also one electron
				if((*iter_mc)->particleProperty()==22&&(*iter_mc)->primaryParticle()){
					for(int ll =0; ll<4;ll++) m_rISR[ll]+=(*iter_mc)->initialFourMomentum()[ll]; 
				}

				//Primary Particle for qqbar, the primaryParticle is True and for other MC, uncertain
				//if (!m_isqqbar){ // for non-qqbar MC
				//	if ((*iter_mc)->primaryParticle()) continue;
				//}
				//if ((*iter_mc)->decayFromGenerator()==0) continue;

				int pid = (*iter_mc)->particleProperty();
				int mother= (*iter_mc)->mother().particleProperty();

				// decayFromGenerator is always true
				//if ((!(*iter_mc)->decayFromGenerator())&&pid!=310 && pid!=130 && mother !=310 &&mother!=130 ) continue;
				if ((*iter_mc)->decayFromGenerator()) n_Gent4++;

				// single electron 
				if (mother==11 && pid==11 && (*iter_mc)->primaryParticle()==1 ) {
					continue;
				}

				//if (abs(pid)==211) cnt_pi++;
				HepLorentzVector  p4_mc_raw = (*iter_mc)->initialFourMomentum();
				m_pid[ntot_MC] = pid;
				m_mother[ntot_MC] = mother;
				m_motherindex[ntot_MC] = (*iter_mc)->mother().trackIndex();
				m_mc_primaryParticle[ntot_MC] = (*iter_mc)->primaryParticle();
				m_mc_leafParticle[ntot_MC] = (*iter_mc)->leafParticle();
				m_mc_decayFromGen[ntot_MC] = (*iter_mc)->decayFromGenerator();
				m_mc_decayInFlight[ntot_MC] = (*iter_mc)->decayInFlight();
				m_mc_children[ntot_MC] = (*iter_mc)->daughterList().size();
				if ((*iter_mc)->daughterList().size()==0) {
					m_mc_childrenFlight[ntot_MC] = 0;
				}else {
					m_mc_childrenFlight[ntot_MC] = (*iter_mc)->daughterList()[0]->decayInFlight();
				}
				m_mc_isFinal[ntot_MC] = false;
				if (m_mc_leafParticle[ntot_MC]==1){
					m_mc_isFinal[ntot_MC] = true;
				} else if (m_mc_leafParticle[ntot_MC]==0 && m_mc_childrenFlight[ntot_MC]==1){
					m_mc_isFinal[ntot_MC] = true;
				} 
				m_pid_2[ntot_MC] = pid;
				m_mother_2[ntot_MC] = mother;
				for(int ll =0; ll<4;ll++) m_p4mc_2[ntot_MC][ll] = p4_mc_raw[ll]; 
				for(int ll =0; ll<4;ll++) m_p4mc[ntot_MC][ll] = p4_mc_raw[ll]; 
				ntot_MC++;

				//Z0
				if(pid==23){/*cout<<"pid is 23"<<endl;*/ 
					for(int ll =0; ll<4;ll++) m_Z0[ll]+=p4_mc_raw[ll]; 
				}

				//91==cluster, 92==string, 94==CMshower
				if(pid==91||pid==92||pid==94){ 
					for(int ll =0; ll<4;ll++){ m_string[ll]=p4_mc_raw[ll];} m_fid = pid;
				}
				//if(pid==22&&(*iter_mc)->primaryParticle()){for(int ll =0; ll<4;ll++) m_rISR[ll]+=p4_mc_raw[ll]; }
				// in qqbar MC mother = pid
				if(mother==91||mother==92||mother==94){ for(int ll =0; ll<4;ll++) m_daughter[ll]+=p4_mc_raw[ll];}

				// no gluo in the qqbar MC
				if(pid==21&&(fabs(mother)==1||fabs(mother)==2||fabs(mother)==3)){
					m_gluon++;
					m_gluon_2++;
					m_gluon_3++;
				}


				// final paricle

				bool permit =false;
				//if(((!(*iter_mc)->primaryParticle())&&(*iter_mc)->leafParticle()&&abs(pid)!=2&&abs(pid)!=3&&abs(pid)!=1&&abs(pid)!=4&&abs(pid)!=5&&abs(pid)!=21)||((*iter_mc)->primaryParticle()&&abs(pid)==22)) permit = true;
				if (m_isqqbar){
					if((!(*iter_mc)->primaryParticle())&&(abs(pid)==11||abs(pid)==13||abs(pid)==211||abs(pid)==321||abs(pid)==2212||abs(pid)==22||abs(pid)==12||abs(pid)==14||abs(pid)==2112))  permit =true;
					if( permit && (abs(mother)==11||abs(mother)==13||abs(mother)==211|abs(mother)==321||abs(mother)==2212||abs(mother)==22))  permit =false;
					if( (abs(pid)==11 || abs(pid)==12 || abs(pid)==13 || abs(pid)==14) &&(abs(mother)==130 ||abs(mother)==321 ||abs(mother)==111 ||abs(mother)==211 ||abs(mother)==13) ){
						permit = false;
					}
					if (abs(pid)==211){
						if (abs(mother)==321) permit=false;
						if (abs(mother)==130) permit=false;
						if (abs(mother)==211) permit=true;
						if (abs(mother)==310) permit=true;
						if (abs(mother)==3122) permit=true;
					}
					if (abs(pid)==2112 && abs(mother==3122)) permit = true; 
					if (abs(pid)==321  && abs(mother==333)) permit = true; 
				} else{
					if(((*iter_mc)->primaryParticle())&&(abs(pid)==11||abs(pid)==13||abs(pid)==211||abs(pid)==321||abs(pid)==2212||abs(pid)==22||abs(pid)==12||abs(pid)==14||abs(pid)==2112))  permit =true;
				}

				if(permit){
					HepLorentzVector  p4_mc= (*iter_mc)->initialFourMomentum();
					if (m_applyBoost) p4_mc.boost(-0.011, 0, 0);
					Hep3Vector p3_mc =  p4_mc.vect();
					m_leafId[nLeaf] = pid;
					m_p4_leaf[nLeaf][0] = p4_mc[0];
					m_p4_leaf[nLeaf][1] = p4_mc[1];
					m_p4_leaf[nLeaf][2] = p4_mc[2];
					m_p4_leaf[nLeaf][3] = p4_mc[3];
					int charge_tmp=-10;
					m_totP4[0]+=p4_mc_raw.px();
					m_totP4[1]+=p4_mc_raw.py();
					m_totP4[2]+=p4_mc_raw.pz();
					m_totP4[3]+=p4_mc_raw.e();
					//if(m_debug) cout<<"using m_particleTable  "<< endl;
					if(pid>0 ) {
						if(m_particleTable->particle( pid )) charge_tmp = m_particleTable->particle( pid )->charge();
					}
					else if(pid<0){
						if(m_particleTable->particle( -pid )) charge_tmp = -1 * m_particleTable->particle( -pid )->charge();
					}
					if (m_debug) cout<<"pid = "<<pid<<" charge = "<<charge_tmp<<endl;
					if(pid==21) charge_tmp = 0;
					if(abs(pid)==22) nLeaf_netural++;
					if(abs(pid)==111) nmc_pi0++;

					//if(m_debug) cout<<"end of using m_particleTable  "<<endl;

					//if(m_debug) cout<<"charge of the mc particle "<< charge_tmp<<endl;

					if(charge_tmp!=0){ //only charge hadrons use to weight 
						mc_charged_hadron.push_back(p4_mc); 
						mc_charge.push_back(charge_tmp);
						mc_charged_trkIndex.push_back((*iter_mc)->trackIndex());
						mc_pid.push_back(pid);
						mc_mother.push_back(mother);
						m_mc_cos[nLeaf_charged] =  p4_mc_raw.vect().cosTheta();
						m_mc_theta[nLeaf_charged] =  p4_mc_raw.vect().theta();
						m_mc_phi[nLeaf_charged] =  p4_mc_raw.vect().phi();
						m_mc_theta_2[nLeaf_charged] =  p4_mc_raw.vect().theta();
						m_mc_cos_2[nLeaf_charged] =  p4_mc_raw.vect().cosTheta();
						m_mc_phi_2[nLeaf_charged] =  p4_mc_raw.vect().phi();
						nLeaf_charged++;
						m_net_leaf_charge_2+= charge_tmp;
						if(abs(pid)==211) cnt_pi++;
					}

					for(int aa=0;aa<Nbin_1;aa++){/// all leaf particles are used to calculate Thrust axis to check
						for(int bb=0;bb<Nbin_2;bb++) {
						double theta_in = 0+step_1*(aa+0.5);
						double phi_in = 0+step_2*(bb+0.5);
						double rx=sin(theta_in)*cos(phi_in);
						double ry=sin(theta_in)*sin(phi_in);
						double rz=-1*cos(theta_in);
						Hep3Vector axis_tmp(rx,ry,rz);
						double ang=axis_tmp.angle(p3_mc);
						Thrust_record_mc[aa][bb]+=fabs(p3_mc.mag()*cos(ang));
						Thrust_unit_mc[aa][bb]+=fabs(p3_mc.mag());
					}
				}
				nLeaf++;
				if(nLeaf>=60) cout<<"warning !!! number of leaf particle exceed 60 !"<<endl;
				if(m_debug) cout<<" number of leaf particle  : "<< nLeaf<<endl;
				}//end of leaf particle


				if(mother==23&&(pid==1|| pid==2||pid==3||pid==4))  { p4_q = p4_mc_raw; cnt_q++;m_qid=pid;m_qid_2=pid;}
				if(mother==23&&(pid==-1|| pid==-2||pid==-3||pid==-4))  { p4_qbar = p4_mc_raw; cnt_qbar++;m_qbarid = pid; m_qbarid_2 = pid;}

				if ((pid)==23) {
					MCmode = m_svc->extract(*iter_mc, pdgid, motherindex);
				}
				m_numParticle += 1;
			}//end of loop of MCParticles
			//if(m_net_leaf_charge_2!=0)  cout<<"something is wrong ! net charge is "<< m_net_leaf_charge_2<<endl;
			// if(nLeaf<1) cout<<" warning: no leaf particle ?!"<<endl;
		}//end of reading MCParticles 

		///booking
		for(int ii=0;ii<4;ii++){ m_p4_q[ii]= p4_q[ii]; m_p4_qbar[ii] = p4_qbar[ii];m_p4_q_2[ii]= p4_q[ii]; m_p4_qbar_2[ii] = p4_qbar[ii]; m_q_p4[ii]=p4_q[ii];m_qbar_p4[ii] = p4_qbar[ii];}
		Thr_axis_truth = p4_q.vect();

		plane_ref_mc=  Thr_axis_truth.cross(epem);
		if(mc_charged_hadron.size()!=mc_charge.size()||mc_charged_hadron.size()!=mc_charged_trkIndex.size())
			if(m_debug) cout<<"warning!  mc_charged_hadron.size() "<< mc_charged_hadron.size() << " , "<< mc_charge.size()<< " , "<< mc_charged_trkIndex.size() <<endl;

		//assign charged hadron-pairs using Thr_axis_truth
		for(int ll=0;ll<mc_charged_hadron.size();ll++) {
			//book
			for (int tt=0;tt<4;tt++)  { m_mom_leaf[ll][tt] =  mc_charged_hadron[ll][tt];  m_mom_leaf_2[ll][tt] =  mc_charged_hadron[ll][tt];}
			m_mom_leaf[ll][4] = mc_charge[ll];
			m_mom_leaf[ll][5] = mc_pid[ll];
			m_mom_leaf[ll][6] = mc_mother[ll];
			m_mom_leaf_2[ll][4] = mc_charge[ll];
			m_mom_leaf_2[ll][5] = mc_pid[ll];
			m_mom_leaf_2[ll][6] = mc_mother[ll];
			m_charge_mc[ll] = mc_charge[ll];
			m_charge_mc_2[ll] = mc_charge[ll];
			m_trkIndex_mc[ll] = mc_charged_trkIndex[ll];

			plane_mc =  mc_charged_hadron[ll].vect().cross(Thr_axis_truth);
			double angle_t=0;
			angle_t= mc_charged_hadron[ll].vect().angle(Thr_axis_truth);
			if(angle_t>-1*(CLHEP::pi)/2&&angle_t<(CLHEP::pi)/2) m_belong_mc[ll]= 1;
			else m_belong_mc[ll]= -1;
			m_belong_mc_2[ll] = m_belong_mc[ll];
			m_z_mc[ll]=  mc_charged_hadron[ll].e()/m_beamE;
			m_z_mc_2[ll] = m_z_mc[ll];
			double an1= mc_charged_hadron[ll].vect().angle(Thr_axis_truth);
			Hep3Vector tmp = (epem.cross(Thr_axis_truth)).cross(Thr_axis_truth.cross(mc_charged_hadron[ll].vect()));
			double sign= Thr_axis_truth.angle(tmp);
			m_pt_mc[ll] = mc_charged_hadron[ll].vect().mag()*cos(an1);
			if(sign<(CLHEP::pi)/2) {
				m_phi_mc[ll] = plane_mc.angle(plane_ref_mc); 
			} else {
				m_phi_mc[ll] = -1*plane_mc.angle(plane_ref_mc);
			}

			for(int tt=ll+1;tt<mc_charged_hadron.size();tt++) {
				m_angle_mc[ll][tt] = mc_charged_hadron[ll].vect().angle(mc_charged_hadron[tt].vect());
				m_mass_mc_2[ll][tt] = (mc_charged_hadron[ll]+mc_charged_hadron[tt]).m();
				m_mass_mc[ll][tt] = (mc_charged_hadron[ll]+mc_charged_hadron[tt]).m();
				m_angle_mc_2[ll][tt] = m_angle_mc[ll][tt];
				m_isfavor_mc[ll][tt] = -1* mc_charge[ll] * mc_charge[tt];
				m_isfavor_mc_2[ll][tt] = m_isfavor_mc[ll][tt];
				///calculate 
				m_qt_mc[ll][tt] =  mc_charged_hadron[ll].vect().mag()*sin(m_angle_mc[ll][tt])/m_z_mc[ll];
				m_qt_mc_2[ll][tt] = m_qt_mc[ll][tt];
				Hep3Vector plane_ref0_mc = epem.cross(mc_charged_hadron[tt].vect());	
				Hep3Vector plane_3_mc   = mc_charged_hadron[tt].vect().cross(mc_charged_hadron[ll].vect());
				Hep3Vector plane_cross = plane_ref0_mc.cross( plane_3_mc);
				double phi0_mc = plane_ref0_mc.angle(plane_3_mc);
				if(mc_charged_hadron[tt].angle(plane_cross)>(CLHEP::pi)/2) phi0_mc *=-1;
				phi0_mc *=2;

				if(phi0_mc<-1*(CLHEP::pi))  phi0_mc+=2*CLHEP::pi;
				else if(phi0_mc>CLHEP::pi) phi0_mc-=2*CLHEP::pi;
				m_phi0_mc [ll][tt] = phi0_mc;
				m_phi0_mc_2[ll][tt] = phi0_mc;
			}
		}
	} 

	for (int i=0; i< pdgid.size(); i++ ) {
		m_pdgid[i] = pdgid[i];
		m_motheridx[i] = motherindex[i];
	}

	m_tot = ntot_MC;
	m_idxmc = pdgid.size();
	m_nleaf  = nLeaf;
	m_nleaf_ch  = nLeaf_charged;
	m_nleaf_ne= nLeaf_netural;
	m_n_pi0= nmc_pi0;

	m_nleaf_ch_2 = nLeaf_charged;
	m_nleaf_ne_2 = nLeaf_netural;
	m_n_pi0_2= nmc_pi0;
	m_npion_mc = cnt_pi;
	m_npion_mc_2 = cnt_pi;
	m_tot_2 = ntot_MC;
	m_nleaf_2  = nLeaf;
	m_n_Gent4 =  n_Gent4;

	if(runNo<0) m_tuple2->write();
	if(runNo<0) m_tuple3->write();

	if(runNo<0) {
		for(int aa=0;aa<Nbin_1;aa++) {
			for(int bb=0;bb<Nbin_2;bb++) {
				double theta_in = 0+step_1*(aa+0.5);
				double phi_in = 0+step_2*(bb+0.5);
				double rx=sin(theta_in)*cos(phi_in);
				double ry=sin(theta_in)*sin(phi_in);
				double rz=-1*cos(theta_in);
				Hep3Vector axis_tmp(rx,ry,rz);
				Thrust_record_mc[aa][bb] = Thrust_record_mc[aa][bb]/Thrust_unit_mc[aa][bb];
				//if(m_debug) cout<<" rx : "<< rx  <<endl;
				//if(m_debug) cout<<" aa: "<<aa << " , bb "<< bb <<"  , "<< Thrust_record_mc[aa][bb]<<endl;
				if(Thrust_record_mc[aa][bb]>Thrust_diy_mc) {
					Thrust_diy_mc=  Thrust_record_mc[aa][bb];
					Thr_axis_mc  = axis_tmp;
				}
			}
		}
	}

	m_Thr_diy_mc = Thrust_diy_mc;
	for(int ii=0;ii<3;ii++) { m_Thr_axis_mc[ii] = Thr_axis_mc[ii];}


	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	if(!evtRecEvent){
		log << MSG::ERROR << "EvtRecEvent not found" << endreq;
		return StatusCode::FAILURE;
	}
	log << MSG::DEBUG <<"ncharg, nneu, tottks = "
		<< evtRecEvent->totalCharged() << " , "
		<< evtRecEvent->totalNeutral() << " , "
		<< evtRecEvent->totalTracks() <<endreq;


	/*---------------------------------------------
					           begin rec
	---------------------------------------------*/
	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);
	int ngood=0; int nshower=0; int npion=0; int nkaon=0;
	Vint iGood;   iGood.clear(); 
	Vint iShower; iShower.clear();
	Vint iPion;   iPion.clear();
	Vint iKaon;   iKaon.clear();
	Vint iEle;    iEle.clear();
	Vint iProton; iProton.clear();
	//Vint iMuon;  iMuon.clear();
	double Evisible=0;

	ISimplePIDSvc* m_simplePIDSvc;
	Gaudi::svcLocator()->service("SimplePIDSvc", m_simplePIDSvc);

	if(m_debug) {
		cout<<"EXE::total charged size : "<< evtRecEvent->totalCharged()<< "  total neutral tracks " << evtRecEvent->totalNeutral() <<endl;
	}

	TrackMatch m_getMatchIndex;
	m_getMatchIndex.setDebug(m_debug);

	if(m_debug){
		cout<<"EXE:: judge good tracks"<<endl;
	}
	for(int i=0; i<evtRecEvent->totalCharged(); i++){
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(Is_good_trk(itTrk)))  continue;
		iGood.push_back(i);
		RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
		RecMdcKalTrack* mdcKalTrk= (*itTrk)->mdcKalTrack();
		HepLorentzVector p4;
		p4=mdcKalTrk->p4(xmass[2]);
		Evisible +=p4.e();
		m_simplePIDSvc->preparePID(*itTrk);
		if(m_simplePIDSvc->iskaon()) {iKaon.push_back(i);}
		if(IsPronton(itTrk))         {iProton.push_back(i);}
		if(m_simplePIDSvc->ispion()) {iPion.push_back(i);}
		if(m_simplePIDSvc->iselectron()) {iEle.push_back(i);}
	}

	if(m_debug){
		cout<<"EXE:: judge good shower"<<endl;
	}
	for(int i=evtRecEvent->totalCharged();i<evtRecEvent->totalTracks();i++){
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!Is_good_gam(itTrk)) continue;
		//nshower++;
		iShower.push_back(i);
	}

	ngood= iGood.size(); 
	nshower=iShower.size();
	npion= iPion.size();
	nkaon = iKaon.size();
	m_npion=npion;
	m_nkaon=nkaon;
	m_nele=iEle.size();
	//m_nmuon=iMuon.size();
	m_nproton=iProton.size();
	m_ngood=ngood;
	//m_ngood_3=ngood;
	m_nshower=nshower;
	if(m_debug) cout<<"npion: "<< npion<<endl;

	if(iEle.size()>m_nele_cut) return StatusCode::SUCCESS;
	if(iGood.size()>=20) {
		if (m_debug) { 
			cout<<"abnormal! ngood == "<< iGood.size()<<endl;
		}
		return StatusCode::SUCCESS; 
	}
	if(ngood<2){
		if(m_debug) {
			cout<<"can't pass ngood "<<ngood<< " < " << 2<<endl; 
		}
		return StatusCode::SUCCESS;
	}
	if(ngood+nshower<numParticle_cut){
		if(m_debug) {
			cout<<"can't pass numParticle_cut "<<npion+nshower<< " < " << numParticle_cut<<endl; 
		}
		return StatusCode::SUCCESS;
	}
	if(npion>20||nshower>50){ cout<<"npion exceed 20"<<endl; return StatusCode::SUCCESS; }
	if(nshower>=30) {
		if (m_debug) { 
			cout<<"abnormal! nshower == "<< nshower <<endl;
		}
		return StatusCode::SUCCESS; 
	}



	Hep3Vector  Thr_axis(-1,-1,-1);
	double Thrust_diy=0;
	double Thrust_record[Nbin_1][Nbin_2];
	double Thrust_unit[Nbin_1][Nbin_2];

	for(int aa=0;aa<Nbin_1;aa++) {
		for(int bb=0;bb<Nbin_2;bb++){ Thrust_record[aa][bb] = 0; Thrust_unit[aa][bb]=0;}
	}


	if (m_debug) cout<<__LINE__<<endl;
	for(int i=0;i<iGood.size();i++){
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + iGood[i];
		HepLorentzVector p4(0,0,0,0);
		RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
		RecMdcKalTrack* mdcKalTrk= (*itTrk)->mdcKalTrack();
		p4=mdcKalTrk->p4(xmass[2]);

		if (m_applyBoost) p4.boost(-0.011, 0, 0);

		for ( int ll = 0; ll < 4; ll++ )  m_fourmom_ch[i][ll] = p4[ll];
		m_mom_ch[i] = p4.vect().mag();///momentum of charged track;
		Hep3Vector p3 =  p4.vect();

		for(int aa=0;aa<Nbin_1;aa++) {
			for(int bb=0;bb<Nbin_2;bb++) {
				double theta_in =  0+step_1*(aa+0.5);
				double phi_in = 0+step_2*(bb+0.5);
				double rx=sin(theta_in)*cos(phi_in);
				double ry=sin(theta_in)*sin(phi_in);
				double rz=-1*cos(theta_in);
				Hep3Vector axis_tmp(rx,ry,rz);
				double ang=axis_tmp.angle(p3);
				Thrust_record[aa][bb]+=fabs(p3.mag()*cos(ang)); 
				Thrust_unit[aa][bb]+=fabs(p3.mag());
			}
		}
	}
	if (m_debug) cout<<__LINE__<<endl;

	for(int i=0; i<iShower.size(); i++){
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + iShower[i];
		RecEmcShower* emcShower = (*itTrk)->emcShower();

		double eraw = emcShower->energy();
		double phi =  emcShower->phi();
		double the =  emcShower->theta();
		HepLorentzVector p4_shower (eraw * sin(the) * cos(phi), eraw * sin(the) * sin(phi), eraw * cos(the), eraw ) ;
		if (m_applyBoost) p4_shower.boost(-0.011, 0, 0);

		for ( int ll = 0; ll < 4; ll++ )  m_fourmom_ne[ngood+i][ll] = p4_shower[ll];
		Hep3Vector p3_shower=p4_shower.vect();

		for(int aa=0;aa<Nbin_1;aa++) {
			for(int bb=0;bb<Nbin_2;bb++) {
				double theta_in =  0+step_1*(aa+0.5);
				double phi_in = 0+step_2*(bb+0.5);
				double rx=sin(theta_in)*cos(phi_in);
				double ry=sin(theta_in)*sin(phi_in);
				double rz=-1*cos(theta_in);
				Hep3Vector axis_tmp(rx,ry,rz);
				double ang=axis_tmp.angle(p3_shower);
				Thrust_record[aa][bb]+=fabs(p3_shower.mag()*cos(ang));
				Thrust_unit[aa][bb] += fabs(p3_shower.mag());
			}
		}
		Evisible += p4_shower.e();
	}
	if (m_debug) cout<<__LINE__<<endl;
	if (Evisible>10) return StatusCode::SUCCESS;
	m_evis=Evisible;
	if (m_debug) cout<<__LINE__<<endl;

	for(int aa=0;aa<Nbin_1;aa++) {
		for(int bb=0;bb<Nbin_2;bb++) {
			double theta_in = 0+step_1*(aa+0.5);
			double phi_in = 0+step_2*(bb+0.5);
			double rx=sin(theta_in)*cos(phi_in);
			double ry=sin(theta_in)*sin(phi_in);
			double rz=-1*cos(theta_in);
			Hep3Vector axis_tmp(rx,ry,rz);
			Thrust_record[aa][bb] = Thrust_record[aa][bb]/Thrust_unit[aa][bb];
			if(Thrust_record[aa][bb]>Thrust_diy) {
				Thrust_diy=  Thrust_record[aa][bb];
				Thr_axis  = axis_tmp;
			}
		}
	}
	if (m_debug) cout<<__LINE__<<endl;

	//book
	m_Thr_diy= Thrust_diy;
	if(m_UseTruthAxis==true)  Thr_axis= Thr_axis_truth;
	for(int ii=0;ii<3;ii++)   m_Thr_axis[ii]=Thr_axis[ii];
	if(m_debug) cout<<"get the Thrust axis: "<< Thr_axis[0] <<" , "<< Thr_axis[1]<< " , "<< Thr_axis[2]<<endl;
	m_Thr_theta= Thr_axis.angle(epem);
	Hep3Vector  plane_ref =  Thr_axis.cross(epem);

	////////////////////PION///////////////////

	Vint charge_pions;     charge_pions.clear();
	Vp4 p4_pions;          p4_pions.clear();
	Vint TrackIndex_match_pions; TrackIndex_match_pions.clear();
	Vint Type_match_pions;       Type_match_pions.clear();

	if (m_debug) cout<<__LINE__<<endl;
	for(int i=0;i<iPion.size();i++){
		if (m_debug)cout<<"line = "<<__LINE__<<endl;
		//prepare the track
		//Fill charge_pions
		//Fill p4_pions, the momentums after boost
		//Fill TrackIndex_match, TrackIndex of TrackMatch
		//Fill Type_match_pions, particle type of TrackMatch
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + iPion[i];
		HepLorentzVector p4(0,0,0,0);
		RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
		RecMdcKalTrack* mdcKalTrk= (*itTrk)->mdcKalTrack();
		p4=mdcKalTrk->p4(xmass[2]);
		if (m_debug)cout<<"line = "<<__LINE__<<endl;

		//m_pion_cosT
		m_pion_cosT[i] = p4.vect().cosTheta();
		//m_pion_Theta
		m_pion_Theta[i] = p4.vect().theta();
		//m_pion_phi
		m_pion_phi[i] = p4.vect().phi();
		if (m_debug)cout<<"line = "<<__LINE__<<endl;

		if (m_applyBoost) p4.boost(-0.011, 0, 0);
		if (m_debug)cout<<"line = "<<__LINE__<<endl;

		//m_pion_belong
		double angle_tmp = p4.vect().angle(Thr_axis);
		//if (m_debug)cout<<"line = "<<__LINE__<<endl;
		//if (m_debug)cout<<"angle_tmp = "<<angle_tmp<<endl;
		//if(angle_tmp>-1*(CLHEP::pi)/2&&angle_tmp<(CLHEP::pi)/2) {
		//	if (m_debug) cout<<"npion1 = "<<m_npion1<<" mpion2 = "<<m_npion2<<endl;
		//	if (m_debug) cout<<"i = "<<i<<endl;
		//	m_pion_belong[i]= 1; 
		//if (m_debug)cout<<"line = "<<__LINE__<<endl;
		//	m_npion1++;
		//	if (m_debug)cout<<"line = "<<__LINE__<<endl;
		//} else {
		//	if (m_debug) cout<<"npion1 = "<<m_npion1<<" mpion2 = "<<m_npion2<<endl;
		//	if (m_debug) cout<<"i = "<<i<<endl;
		//	m_pion_belong[i]= -1;
		//if (m_debug)cout<<"line = "<<__LINE__<<endl;
		//	m_npion2++;
		//	if (m_debug)cout<<"line = "<<__LINE__<<endl;
		//}

		if (m_debug)cout<<"line = "<<__LINE__<<endl;
		//m_pion_charge
		charge_pions.push_back(mdcKalTrk->charge());
		m_pion_charge[i] =  charge_pions[i];

		//m_pion_theta0
		m_pion_theta0[i] = p4.vect().angle(epem);

		//m_pion_pt
		m_pion_pt[i] = p4.vect().perp();

		//type match
		//m_pion_trkIndex
		//m_pion_type_match
		if (m_debug)cout<<"line = "<<__LINE__<<endl;
		m_getMatchIndex.IfType(false);
		int trk_index =  -200;
		if(m_debug) cout<<" begin track matching "<<endl;
		int pitype=2;
		m_getMatchIndex.initialize();
		if(runNo<0) {trk_index=m_getMatchIndex.GetTrackIndex(itTrk,pitype);} 
		TrackIndex_match_pions.push_back(trk_index);
		Type_match_pions.push_back(m_getMatchIndex.GetType());// electron 0, muon 1, pion 2, kaon 3, proton 4;
		if(m_debug) cout<<"match_ed   trk_index: "<< trk_index<<endl;
		m_pion_trkIndex[i] =  TrackIndex_match_pions[i];
		m_pion_type_match[i] =  Type_match_pions[i];

		//m_pion_z
		m_pion_z[i] = p4.e()/m_beamE;

		//m_pion_phi1
		Hep3Vector tmp = (epem.cross(Thr_axis)).cross(Thr_axis.cross(p4.vect()));
		Hep3Vector plane = p4.vect().cross(Thr_axis);
		double sign= Thr_axis.angle(tmp);
		if(sign<(CLHEP::pi/2)) {
			m_pion_phi1[i] = plane.angle(plane_ref);
		} else {
			m_pion_phi1[i] = -1*plane.angle(plane_ref); 
		}

		//m_pion_p4;
		p4_pions.push_back(p4);
		for (int jj = 0; jj<4; jj++) m_pion_p4[i][jj] = p4[jj];
	}

	if(TrackIndex_match_pions.size()!=npion)cout<<"waring !"<<endl;;

	///assige pions into hemisphere

	int  com_cnt=0;

	for(int i=0;i<p4_pions.size();i++) {
		HepLorentzVector p4_i = p4_pions[i];
		for(int j=i+1;j<p4_pions.size();j++){
			HepLorentzVector p4_j = p4_pions[j];

			//m_pion_angle
			double angle_tmp=p4_i.vect().angle(p4_j.vect());
			m_pion_angle[i][j] = angle_tmp;

			//m_pion_mass
			double mass = m_pion_mass[i][j]=(p4_i+p4_j).m();
			m_pion_mass[i][j] = mass;

			//m_pion_isfavor
			m_pion_isfavor[i][j] = -1* charge_pions[i] * charge_pions[j];
			//m_pion_ptvs
			m_pion_ptvs[i][j] =   p4_i.vect().mag()*sin(angle_tmp);
			//m_pion_qtvs
			m_pion_qtvs[i][j] =  m_pion_ptvs[i][j]/m_pion_z[i];
			//m_pion_phi0
			m_pion_phi0[i][j] = GetPhi0(p4_i, p4_j);
			if (m_debug){
				cout<<"pioni = "<<i<<"\t"<<"p4 = "<<p4_i.x()<<"\t"<<p4_i.y()<<"\t"<<p4_i.z()<<"\t"<<p4_i.e()<<endl;
				cout<<"pionj = "<<j<<"\t"<<"p4 = "<<p4_j.x()<<"\t"<<p4_j.y()<<"\t"<<p4_j.z()<<"\t"<<p4_j.e()<<endl;
				cout<<"phi0 = "<<m_pion_phi0[i][j]<<endl;
			}
			com_cnt++;
		}
	}

	if(m_debug) cout<<" number of pion combination: "<< com_cnt<<endl;

	////////////////////KAON///////////////////

	Vint charge_kaons;     charge_kaons.clear();
	Vp4 p4_kaons;          p4_kaons.clear();
	Vint TrackIndex_match_kaons; TrackIndex_match_kaons.clear();
	Vint Type_match_kaons;       Type_match_kaons.clear();

	for(int i=0;i<iKaon.size();i++){
		//prepare the track
		//Fill charge_kaons
		//Fill p4_kaons, the momentums after boost
		//Fill TrackIndex_match_kaons, TrackIndex of TrackMatch
		//Fill Type_match_kaons, particle type of TrackMatch
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + iKaon[i];
		HepLorentzVector p4(0,0,0,0);
		RecMdcKalTrack::setPidType(RecMdcKalTrack::kaon);
		RecMdcKalTrack* mdcKalTrk= (*itTrk)->mdcKalTrack();
		p4=mdcKalTrk->p4(xmass[3]);

		//m_kaon_cosT
		m_kaon_cosT[i] = p4.vect().cosTheta();
		//m_kaon_Theta
		m_kaon_Theta[i] = p4.vect().theta();
		//m_kaon_phi
		m_kaon_phi[i] = p4.vect().phi();

		if (m_applyBoost) p4.boost(-0.011, 0, 0);

		//m_kaon_belong
		double angle_tmp = p4.vect().angle(Thr_axis);
		if(angle_tmp>-1*(CLHEP::pi)/2&&angle_tmp<(CLHEP::pi)/2) {
			m_kaon_belong[i]= 1; 
			m_nkaon1++;
		} else {
			m_kaon_belong[i]= -1;
			m_nkaon2++;
		}

		//m_kaon_charge
		charge_kaons.push_back(mdcKalTrk->charge());
		m_kaon_charge[i] =  charge_kaons[i];

		//m_kaon_theta0
		m_kaon_theta0[i] = p4.vect().angle(epem);

		//m_kaon_pt
		m_kaon_pt[i] = p4.vect().perp();

		//type match
		//m_kaon_trkIndex
		//m_kaon_type_match
		m_getMatchIndex.IfType(false);
		int trk_index =  -200;
		if(m_debug) cout<<" begin track matching "<<endl;
		int pitype=3;
		m_getMatchIndex.initialize();
		if(runNo<0) {trk_index=m_getMatchIndex.GetTrackIndex(itTrk,pitype);} 
		TrackIndex_match_kaons.push_back(trk_index);
		Type_match_kaons.push_back(m_getMatchIndex.GetType());// electron 0, muon 1, kaon 2, kaon 3, proton 4;
		if(m_debug) cout<<"match_ed   trk_index: "<< trk_index<<endl;
		m_kaon_trkIndex[i] =  TrackIndex_match_kaons[i];
		m_kaon_type_match[i] =  Type_match_kaons[i];

		//m_kaon_z
		m_kaon_z[i] = p4.e()/m_beamE;

		//m_kaon_phi1
		Hep3Vector tmp = (epem.cross(Thr_axis)).cross(Thr_axis.cross(p4.vect()));
		Hep3Vector plane = p4.vect().cross(Thr_axis);
		double sign= Thr_axis.angle(tmp);
		if(sign<(CLHEP::pi/2)) {
			m_kaon_phi1[i] = plane.angle(plane_ref);
		} else {
			m_kaon_phi1[i] = -1*plane.angle(plane_ref); 
		}

		//m_kaon_p4;
		p4_kaons.push_back(p4);
		for (int jj = 0; jj<4; jj++) m_kaon_p4[i][jj] = p4[jj];
	}

	if(TrackIndex_match_kaons.size()!=nkaon)cout<<"waring !"<<endl;;

	///assige kaons into hemisphere

	com_cnt=0;

	for(int i=0;i<p4_kaons.size();i++) {
		HepLorentzVector p4_i = p4_kaons[i];
		for(int j=i+1;j<p4_kaons.size();j++){
			HepLorentzVector p4_j = p4_kaons[j];

			//m_kaon_angle
			double angle_tmp=p4_i.vect().angle(p4_j.vect());
			m_kaon_angle[i][j] = angle_tmp;

			//m_kaon_mass
			double mass = m_kaon_mass[i][j]=(p4_i+p4_j).m();
			m_kaon_mass[i][j] = mass;

			//m_kaon_isfavor
			m_kaon_isfavor[i][j] = -1* charge_kaons[i] * charge_kaons[j];
			//m_kaon_ptvs
			m_kaon_ptvs[i][j] =   p4_i.vect().mag()*sin(angle_tmp);
			//m_kaon_qtvs
			m_kaon_qtvs[i][j] =  m_kaon_ptvs[i][j]/m_kaon_z[i];

			m_kaon_phi0[i][j] = GetPhi0(p4_i, p4_j);
			if (m_debug){
				cout<<"kaoni = "<<i<<"\t"<<"p4 = "<<p4_i.x()<<"\t"<<p4_i.y()<<"\t"<<p4_i.z()<<"\t"<<p4_i.e()<<endl;
				cout<<"kaonj = "<<j<<"\t"<<"p4 = "<<p4_j.x()<<"\t"<<p4_j.y()<<"\t"<<p4_j.z()<<"\t"<<p4_j.e()<<endl;
				cout<<"phi0 = "<<m_kaon_phi0[i][j]<<endl;
			}
			com_cnt++;
		}
	}

	if(m_debug) cout<<" number of kaon combination: "<< com_cnt<<endl;


	////////////////////KPI///////////////////
	
	for (int i = 0; i< iPion.size(); i++){
		HepLorentzVector p4_i = p4_pions[i];
		for (int j = 0; j< iKaon.size(); j++){
			HepLorentzVector p4_j = p4_kaons[j];
			
			//m_kpi_angle
			double angle_tmp=p4_i.vect().angle(p4_j.vect());
			m_kpi_angle[j][i] = angle_tmp;
			//m_kpi_mass
			double mass = m_kpi_mass[j][i]=(p4_i+p4_j).m();
			m_kpi_mass[j][i] = mass;
			//m_kpi_isfavor
			m_kpi_isfavor[j][i] = -1* charge_kaons[j] * charge_pions[i];
			//m_kpi_ptvs
			m_kpi_ptvs[j][i] =   p4_j.vect().mag()*sin(angle_tmp);
			//m_kpi_qtvs
			m_kpi_qtvs[j][i] =  m_kpi_ptvs[j][i]/m_kaon_z[j];
			//m_pik_ptvs
			m_pik_ptvs[i][j] =   p4_i.vect().mag()*sin(angle_tmp);
			//m_pik_qtvs
			m_pik_qtvs[j][i] =  m_pik_ptvs[i][j]/m_pion_z[i];

			m_kpi_phi0[j][i] = GetPhi0(p4_j, p4_i);
			m_pik_phi0[i][j] = GetPhi0(p4_i, p4_j);

		}
	}

	m_cout_pass++;
	if (m_debug)cout<<"*****************m_cout_pass = "<<m_cout_pass<<endl;
	///retrive MC Truth info

	for(int ll =0 ;ll<mc_charged_trkIndex.size();ll++) {
		bool rec_succ=false;
		for(int tt =0 ; tt<npion;tt++) {
			if (TrackIndex_match_pions[tt] == mc_charged_trkIndex[ll] )  rec_succ =true;
		}
		if(rec_succ) m_Ifrec[ll] = 1;
		else m_Ifrec[ll] = -1;
	}

	m_tuple1->write(); 

	return StatusCode::SUCCESS;
}

StatusCode CollAnaAlg::endRun(){
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"CollAnaAlg::endRun()"<<endreq;
	//add your code here
	return StatusCode::SUCCESS;
}

StatusCode CollAnaAlg::finalize(){
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"CollAnaAlg::finalize()"<<endreq;
	cout<<"CollAnaAlg::total events: "<<  m_cout_all<<endl;
	cout<<"CollAnaAlg::passed events: "<< m_cout_pass<<endl;
	//add your code here
	return StatusCode::SUCCESS;
}


bool CollAnaAlg::Is_good_trk(EvtRecTrackIterator itTrk) {
	double CosThetaCut= m_trk_cos_cut;
	if ( !(*itTrk)->isMdcTrackValid() ) return false;
	if ( !(*itTrk)->isMdcKalTrackValid() ) return false;
	RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();

	Hep3Vector xorigin(0,0,0);
	IVertexDbSvc*  vtxsvc;
	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
	if(vtxsvc->isVertexValid()){
		double* dbv = vtxsvc->PrimaryVertex();
		double* vv = vtxsvc->SigmaPrimaryVertex();
		//pretect
		if(vtxsvc->PrimaryVertex()[0]>100||vtxsvc->PrimaryVertex()[1]>100 || vtxsvc->PrimaryVertex()[2]>100) {
			cout<<"Vertex is abnormal! check your jobOption "<<endl;
			dbv[0]=0;
			dbv[1]=0;
			dbv[2]=0;
		}
		if(m_debug)cout<<dbv[0]<< " , "<< dbv[1] <<" , "<< dbv[2] <<endl;
		xorigin.setX(dbv[0]);
		xorigin.setY(dbv[1]);
		xorigin.setZ(dbv[2]);
	}
	else { cout<<"warning !!! IVertexDbSvc is inValid!"<<endl;}
	if(m_debug) cout<<" xorigin : "<< xorigin[0] << " , "<< xorigin[1]<< " , "<< xorigin[2]<<endl;

	HepVector a = mdcTrk->helix();
	HepSymMatrix Ea = mdcTrk->err();
	HepPoint3D point0(0.,0.,0.);
	HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
	VFHelix helixip3(point0,a,Ea);
	helixip3.pivot(IP);
	HepVector  vecipa = helixip3.a();
	double dr=fabs(vecipa[0]);
	double dz=fabs(vecipa[3]);
	double costheta=cos(mdcTrk->theta());
	if (m_debug)cout<<"dr = "<<dr<<"\t"<<"dz = "<<dz<<endl;
	if (  dr>= m_vr0cut) return false;
	if (  dz>= m_vz0cut ) return false;
	if ( fabs(costheta) >= CosThetaCut ) return false;
	return true;
}

bool CollAnaAlg::Is_good_gam(EvtRecTrackIterator itTrk){
	SmartDataPtr<EvtRecEvent> recEvt(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	SmartDataPtr<EvtRecTrackCol> recTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);
	if(!(*itTrk)->isEmcShowerValid()) return false;

	RecEmcShower *emcTrk = (*itTrk)->emcShower();
	// double  m_minEnergy = 0.025;
	bool  m_useBarrelEndcap   = true;
	double  m_maxCosThetaBarrel = 0.8;
	double  m_minCosThetaEndcap = 0.84;
	double  m_maxCosThetaEndcap = 0.92;
	double  m_minTime      = 0.;
	double  m_maxTime      = 14.;

	double eraw = emcTrk->energy();
	double phi =  emcTrk->phi();
	double the =  emcTrk->theta();
	HepLorentzVector shP4( eraw * sin(the) * cos(phi),
			eraw * sin(the) * sin(phi),
			eraw * cos(the),
			eraw );
	double cosThetaSh = shP4.vect().cosTheta();
	/// Minimum energy
	if (shP4.e() <= m_energy_cut_b) return( false );
	/// Barrel/Endcap
	bool inBarrelEndcap = false;
	if(fabs(cosThetaSh) < m_maxCosThetaBarrel) inBarrelEndcap = true;
	if((fabs(cosThetaSh) > m_minCosThetaEndcap)
			&&(fabs(cosThetaSh) < m_maxCosThetaEndcap)
			&&(shP4.e() > m_energy_cut_e )) inBarrelEndcap = true;
	if(m_useBarrelEndcap&&!inBarrelEndcap) return( false );
	/// Time, only apply timing cuts if "recEvt->totalCharged() > 0"
	if ( m_applyTimeCut ) {
		double time = emcTrk->time();
		if ( recEvt->totalCharged() > 0 ) {
			if ( time < m_minTime || time > m_maxTime ) return false;
		}
		else {
			RecEmcShower* firstG = (*(recTrkCol->begin()))->emcShower();
			double deltaTime = fabs(time - firstG->time());
			if ( deltaTime > 10 ) return false;

		}

	}

	/// Dang
	Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
	if (recEvt->totalCharged() > 0) {
		double dang = 200.;
		for (int j = 0; j < recEvt->totalCharged(); j++) {
			EvtRecTrackIterator jtTrk = recTrkCol->begin() + j;
			if ( !(*jtTrk)->isExtTrackValid() ) continue;
			RecExtTrack* extTrk = (*jtTrk)->extTrack();
			if ( extTrk->emcVolumeNumber() == -1 ) continue;
			Hep3Vector extpos = extTrk->emcPosition();
			double  angd1 = extpos.angle(emcpos);
			if ( angd1 < dang ) dang = angd1;
		}
		if ( dang < 200 ) {
			dang = dang * 180 / (CLHEP::pi);
			if (m_applyDangCut&&(dang <= m_dang_cut)) return( false );
		}
	}  // End of "recEvt->totalCharged() > 0" IF
	return( true );
}

bool CollAnaAlg::IsPronton(EvtRecTrackIterator itTrk) {
	double m_prob_cut =0.001;
	if(!(*itTrk)->isMdcKalTrackValid()) return false;
	ParticleID *pid = ParticleID::instance();
	pid->init();
	pid->setMethod(pid->methodProbability());
	pid->setChiMinCut(4);
	pid->setRecTrack(*itTrk);
	pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2());
	pid->identify(pid->onlyProton()|pid->onlyPion()|pid->onlyKaon());
	pid->calculate();
	if(!(pid->IsPidInfoValid())) return false;

	if(pid->prob(4)<m_prob_cut) return  false;

	if( pid->prob(4)< pid->prob(3) || pid->prob(4)< pid->prob(2)) return false;
	return true;
}

double CollAnaAlg::GetPhi0(HepLorentzVector p4_i, HepLorentzVector p4_j){
	Hep3Vector epem(0,0,1);
	Hep3Vector 	plane_ref_0_b = epem.cross(p4_j.vect());
	Hep3Vector plane_3_b = p4_j.vect().cross(p4_i.vect());
	double tmp_phi0 = plane_3_b.angle(plane_ref_0_b);
	Hep3Vector tmp = plane_ref_0_b.cross(plane_3_b);
	double sign = p4_j.angle(tmp) ;
	if(sign>(CLHEP::pi/2)) tmp_phi0*=-1;
	tmp_phi0*=2;
	if(tmp_phi0<-1*CLHEP::pi)  tmp_phi0 =  tmp_phi0 + 2*CLHEP::pi;
	else if (tmp_phi0>CLHEP::pi) tmp_phi0 = tmp_phi0 - 2*CLHEP::pi;
	return tmp_phi0;
}
