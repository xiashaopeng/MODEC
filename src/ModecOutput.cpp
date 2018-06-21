#include "ModecClass.h"

const double Avogadro_Constant = 6.022140857E+23; // 阿伏加德罗常数   Source: 2014 CODATA
const double Electron_Coulomb = 1.6021766208E-19;//1.6021766208E-19; // 电子电量        Source: 2014 CODATA

const double cutoff = 1e-50; // 浓度的截断误差

using namespace std;
using namespace tinyxml2;

void ModecClass::ModecOutputXml() {
    int nucl_size = ModecNuclideLibrary.nuclide_number_;

    XMLDocument doc;
    doc.Parse("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
    XMLElement* root = doc.NewElement("MODEC");
    doc.InsertEndChild(root);
    root = doc.RootElement();

    int size = evolution_mode_.size();
    int size_tot = 0;
    vector<string> time_unit;
    for(int i = 0; i < size; ++i) {
        size_tot += substep_[i];
    }

    vector<double> cum_time;
    cum_time.resize(size_tot + 1);

    int count = 0;
    for (int i=0; i < size; ++i) {
        for (int j = 0; j < substep_[i]; ++j) {
            cum_time[count + 1] = cum_time[count] + burnup_time_[i];
            time_unit.push_back(time_unit_[i]);
            count++;
        }
    }

    if( if_flow_mode_ == 0 ) {
        XMLElement *time = doc.NewElement("Time");
        time->SetAttribute("unit","s");
        if(print_mode_ == 0) {
            stringstream time_string;
            time_string << scientific << setprecision(output_precision_);
            time_string << cum_time[0];
            time_string << " ";
            time_string << cum_time.back();

            string time_s(time_string.str());
            XMLText *time_text = doc.NewText(time_s.c_str());
            time->InsertEndChild(time_text);
            root->InsertEndChild(time);
        } else {
            stringstream time_string;
            time_string << scientific << setprecision(output_precision_);
            for (int j = 0; j < size_tot; ++j) {
                time_string << cum_time[j] << " ";
            }
            time_string << cum_time.back();

            string time_s(time_string.str());
            //time_string >> time_s;
            XMLText *time_text = doc.NewText(time_s.c_str());
            time->InsertEndChild(time_text);
            root->InsertEndChild(time);
        }

        XMLElement *power = doc.NewElement("Power");
        power->SetAttribute("unit","MW");
        if(print_mode_ == 0) {
            stringstream power_string;
            power_string << scientific << setprecision(output_precision_);
            power_string << power_vector_[0];
            power_string << " ";
            power_string << power_vector_.back();

            string power_s(power_string.str());
            XMLText * power_text = doc.NewText(power_s.c_str());
            power->InsertEndChild(power_text);
            root->InsertEndChild(power);
        } else {
            stringstream power_string;
            power_string << scientific << setprecision(output_precision_);
            for (int j = 0; j < size_tot; ++j) {
                power_string << power_vector_[j] << " ";
            }
            power_string << power_vector_.back();

            string power_s(power_string.str());
            XMLText *power_text = doc.NewText(power_s.c_str());
            power->InsertEndChild(power_text);
            root->InsertEndChild(power);
        }

        XMLElement *flux = doc.NewElement("Flux");
        flux->SetAttribute("unit","cm-2*s-1");
        if(print_mode_ == 0) {
            stringstream flux_string;
            flux_string << scientific << setprecision(output_precision_);
            flux_string << flux_vector_[0];
            flux_string << " ";
            flux_string << flux_vector_.back();

            string flux_s(flux_string.str());
            XMLText * flux_text = doc.NewText(flux_s.c_str());
            flux->InsertEndChild(flux_text);
            root->InsertEndChild(flux);
        } else {
            stringstream flux_string;
            flux_string << scientific << setprecision(output_precision_);
            for (int j = 0; j < size_tot; ++j) {
                flux_string << flux_vector_[j] << " ";
            }
            flux_string << flux_vector_.back();

            string flux_s(flux_string.str());
            XMLText *flux_text = doc.NewText(flux_s.c_str());
            flux->InsertEndChild(flux_text);
            root->InsertEndChild(flux);
        }

        vector<double> kinf;
        vector<double> prod_neu;
        vector<double> abs_neu;
        if (if_print_kinf_ != 0 || if_print_fission_rate_ != 0 || if_print_absorption_rate_ != 0) {
            kinf.resize(size_tot + 1);
            prod_neu.resize(size_tot + 1);
            abs_neu.resize(size_tot + 1);

            for (int j = 0; j <= size_tot; ++j) {
                for (int i = 0; i < nucl_size; ++i) {
                    prod_neu[j] += n_vector_[j][i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24;
                    abs_neu[j] += n_vector_[j][i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24;
                }

                prod_neu[j] *= flux_vector_[j];
                abs_neu[j] *= flux_vector_[j];

                if (abs_neu[j] != 0.0) {
                    kinf[j] = prod_neu[j] / abs_neu[j];
                }
            }

            if (if_print_kinf_ == true) {
                XMLElement* kinf_node = doc.NewElement("K-infinite");
                if (print_mode_ == 0) {
                    stringstream kinf_string;
                    kinf_string << scientific << setprecision(output_precision_);
                    kinf_string << kinf[0];
                    kinf_string << " ";
                    kinf_string << kinf.back();

                    string kinf_s(kinf_string.str());

                    XMLText *kinf_text = doc.NewText(kinf_s.c_str());
                    kinf_node->InsertEndChild(kinf_text);
                    root->InsertEndChild(kinf_node);
                } else {
                    stringstream kinf_string;
                    kinf_string << scientific << setprecision(output_precision_);
                    for (int j = 0; j < size_tot; ++j) {
                        kinf_string << kinf[j] << " ";
                    }
                    kinf_string << kinf.back();

                    string kinf_s(kinf_string.str());
                    XMLText *kinf_text = doc.NewText(kinf_s.c_str());
                    kinf_node->InsertEndChild(kinf_text);
                    root->InsertEndChild(kinf_node);
                }
            }

        }

        XMLElement *nuclides = doc.NewElement("Nuclides");
        nuclides->SetAttribute("zone","core");
        root->InsertEndChild(nuclides);

        XMLElement *concentration = doc.NewElement("Concentration");
        concentration->SetAttribute("unit",dens_unit_.c_str());
        nuclides->InsertEndChild(concentration);

        double convert_coeff;

        for ( int i = 0; i < nucl_size; ++i) {
            XMLElement *nuclide = doc.NewElement("nuclide");
            int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
            int NZ = nucl_id / 10000;
            int NA = (nucl_id - NZ * 10000) / 10;
            int NG = nucl_id - nucl_id / 10 * 10;
            string nucl_name;
            nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
            if(NG == 1) {
                nucl_name += "m";
            }

            if (dens_unit_ == "mol") {
                convert_coeff = 1.0;
            } else if (dens_unit_ == "g") {
                convert_coeff = double(NA);
            } else if (dens_unit_ == "kg") {
                convert_coeff = double(NA) / 1000.0;
            } else if (dens_unit_ == "atom") {
                convert_coeff = Avogadro_Constant;
            } else if (dens_unit_ == "atom/(barn-cm)") {
                convert_coeff = Avogadro_Constant * 1.0e-24;
            }

            nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
            nuclide->SetAttribute("name",nucl_name.c_str());
            if (print_mode_ == 0) {
                stringstream nuclide_string;
                nuclide_string << scientific << setprecision(output_precision_);
                nuclide_string << n_vector_[0][i] * convert_coeff;
                nuclide_string << " ";
                if(n_vector_[size_tot][i] < cutoff) {
                    nuclide_string << 0.0;
                } else {
                    nuclide_string << n_vector_[size_tot][i] * convert_coeff;
                }
                string nuclide_s(nuclide_string.str());
                XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                nuclide->InsertEndChild(nuclide_text);
                concentration->InsertEndChild(nuclide);
            } else {
                stringstream nuclide_string;
                nuclide_string << scientific << setprecision(output_precision_);
                for (int j = 0; j < size_tot; ++j) {
                    if (n_vector_[j][i] < cutoff) {
                        nuclide_string << 0.0 << " ";
                    } else {
                        nuclide_string << n_vector_[j][i] * convert_coeff << " ";
                    }
                }
                if (n_vector_[size_tot][i] < cutoff) {
                    nuclide_string << 0.0 ;
                } else {

                    nuclide_string << n_vector_[size_tot][i] * convert_coeff;
                }
                string nuclide_s(nuclide_string.str());
                XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                nuclide->InsertEndChild(nuclide_text);
                concentration->InsertEndChild(nuclide);
            }

        }
        if (if_print_activity_ == true) {
            XMLElement * act_node = doc.NewElement("RadioActivity");
            act_node->SetAttribute("unit","Ci");
            nuclides->InsertEndChild(act_node);

            vector< double > lambda;
            lambda = ModecNuclideLibrary.nuclide_library_vector_[1];

            for (int i = 0; i < nucl_size; ++ i) {
                XMLElement *nuclide = doc.NewElement("nuclide");
                int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                int NZ = nucl_id / 10000;
                int NA = (nucl_id - NZ * 10000) / 10;
                int NG = nucl_id - nucl_id / 10 * 10;
                string nucl_name;
                nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
                if(NG == 1) {
                    nucl_name += "m";
                }

                nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
                nuclide->SetAttribute("name",nucl_name.c_str());
                if(print_mode_ == 0) {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    nuclide_string << n_vector_[0][i] * Avogadro_Constant * lambda[i] / (3.7e10);
                    nuclide_string << " ";
                    if (n_vector_[size][i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][i] * Avogadro_Constant * lambda[i] / (3.7e10);
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    act_node->InsertEndChild(nuclide);
                } else {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    for ( int j = 0; j < size_tot; ++j) {
                        if(n_vector_[j][i] < cutoff) {
                            nuclide_string << 0.0 <<" ";
                        } else {
                            nuclide_string << n_vector_[j][i] * Avogadro_Constant * lambda[i] / (3.7e10) << " ";
                        }
                    }
                    if (n_vector_[size_tot][i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][i] * Avogadro_Constant * lambda[i] / (3.7e10);
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    act_node->InsertEndChild(nuclide);
                }
            }

        }

        if (if_print_decayenergy_ == true) {

            XMLElement * decay_node = doc.NewElement("DecayHeat");
            decay_node->SetAttribute("unit","Watts");
            nuclides->InsertEndChild(decay_node);

            vector< double > lambda;
            lambda = ModecNuclideLibrary.nuclide_library_vector_[1];

            vector< double > Q;
            Q = ModecNuclideLibrary.nuclide_library_vector_[3];

            for (int i = 0; i < nucl_size; ++ i) {
                XMLElement *nuclide = doc.NewElement("nuclide");
                int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                int NZ = nucl_id / 10000;
                int NA = (nucl_id - NZ * 10000) / 10;
                int NG = nucl_id - nucl_id / 10 * 10;
                string nucl_name;
                nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
                if(NG == 1) {
                    nucl_name += "m";
                }

                nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
                nuclide->SetAttribute("name",nucl_name.c_str());
                if(print_mode_ == 0) {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    nuclide_string << n_vector_[0][i] * Avogadro_Constant * lambda[i] * Q[i] * Electron_Coulomb * 1.0e6;
                    nuclide_string << " ";
                    if (n_vector_[size][i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][i] * Avogadro_Constant * lambda[i] * Q[i] * Electron_Coulomb * 1.0e6;
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    decay_node->InsertEndChild(nuclide);
                } else {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    for ( int j = 0; j < size_tot; ++j) {
                        if(n_vector_[j][i] < cutoff) {
                            nuclide_string << 0.0 <<" ";
                        } else {
                            nuclide_string << n_vector_[j][i] * Avogadro_Constant * lambda[i] * Q[i] * Electron_Coulomb * 1.0e6 << " ";
                        }
                    }
                    if (n_vector_[size_tot][i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][i] * Avogadro_Constant * lambda[i] * Q[i] * Electron_Coulomb * 1.0e6;
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    decay_node->InsertEndChild(nuclide);
                }
            }

        }

        if (if_print_ampc_ == true) {

            XMLElement * ampc_node = doc.NewElement("AMPC");
            ampc_node->SetAttribute("unit","m3-air/Bq");
            nuclides->InsertEndChild(ampc_node);

            vector< double > lambda;
            lambda = ModecNuclideLibrary.nuclide_library_vector_[1];

            vector< double > AMPC;
            AMPC = ModecNuclideLibrary.nuclide_library_vector_[4];

            for (int i = 0; i < nucl_size; ++ i) {
                XMLElement *nuclide = doc.NewElement("nuclide");
                int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                int NZ = nucl_id / 10000;
                int NA = (nucl_id - NZ * 10000) / 10;
                int NG = nucl_id - nucl_id / 10 * 10;
                string nucl_name;
                nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
                if(NG == 1) {
                    nucl_name += "m";
                }

                nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
                nuclide->SetAttribute("name",nucl_name.c_str());
                if(print_mode_ == 0) {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);
                    nuclide_string << n_vector_[0][i] * Avogadro_Constant * lambda[i] / 3.7e10 / AMPC[i];
                    nuclide_string << " ";
                    if (n_vector_[size][i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][i] * Avogadro_Constant * lambda[i] / 3.7e10 / AMPC[i];
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    ampc_node->InsertEndChild(nuclide);
                } else {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    for ( int j = 0; j < size_tot; ++j) {
                        if(n_vector_[j][i] < cutoff) {
                            nuclide_string << 0.0 <<" ";
                        } else {
                            nuclide_string << n_vector_[j][i] * Avogadro_Constant * lambda[i] / 3.7e10 / AMPC[i] << " ";
                        }
                    }
                    if (n_vector_[size_tot][i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][i] * Avogadro_Constant * lambda[i] / 3.7e10 / AMPC[i];
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    ampc_node->InsertEndChild(nuclide);
                }
            }

        }

        if (if_print_wmpc_ == true) {

            XMLElement * wmpc_node = doc.NewElement("WMPC");
            wmpc_node->SetAttribute("unit","m3-water/Bq");
            nuclides->InsertEndChild(wmpc_node);

            vector< double > lambda;
            lambda = ModecNuclideLibrary.nuclide_library_vector_[1];

            vector< double > WMPC;
            WMPC = ModecNuclideLibrary.nuclide_library_vector_[5];

            for (int i = 0; i < nucl_size; ++ i) {
                XMLElement *nuclide = doc.NewElement("nuclide");
                int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                int NZ = nucl_id / 10000;
                int NA = (nucl_id - NZ * 10000) / 10;
                int NG = nucl_id - nucl_id / 10 * 10;
                string nucl_name;
                nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
                if(NG == 1) {
                    nucl_name += "m";
                }

                nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
                nuclide->SetAttribute("name",nucl_name.c_str());
                if(print_mode_ == 0) {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    nuclide_string << n_vector_[0][i] * Avogadro_Constant * lambda[i] / 3.7e10 / WMPC[i];
                    nuclide_string << " ";
                    if (n_vector_[size][i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][i] * Avogadro_Constant * lambda[i] / 3.7e10 / WMPC[i];
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    wmpc_node->InsertEndChild(nuclide);
                } else {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    for ( int j = 0; j < size_tot; ++j) {
                        if(n_vector_[j][i] < cutoff) {
                            nuclide_string << 0.0 <<" ";
                        } else {
                            nuclide_string << n_vector_[j][i] * Avogadro_Constant * lambda[i] / 3.7e10 / WMPC[i] << " ";
                        }
                    }
                    if (n_vector_[size_tot][i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][i] * Avogadro_Constant * lambda[i] / 3.7e10 / WMPC[i];
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    wmpc_node->InsertEndChild(nuclide);
                }
            }

        }

        if (if_print_toxicity_ == true) {

            XMLElement * toxicity_node = doc.NewElement("RadioToxicity");
            toxicity_node->SetAttribute("unit","Sv");
            nuclides->InsertEndChild(toxicity_node);

            vector< double > lambda;
            lambda = ModecNuclideLibrary.nuclide_library_vector_[1];

            vector< double > TOXCI;
            TOXCI = ModecNuclideLibrary.nuclide_library_vector_[10];

            for (int i = 0; i < nucl_size; ++ i) {
                XMLElement *nuclide = doc.NewElement("nuclide");
                int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                int NZ = nucl_id / 10000;
                int NA = (nucl_id - NZ * 10000) / 10;
                int NG = nucl_id - nucl_id / 10 * 10;
                string nucl_name;
                nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
                if(NG == 1) {
                    nucl_name += "m";
                }

                nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
                nuclide->SetAttribute("name",nucl_name.c_str());
                if(print_mode_ == 0) {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    nuclide_string << n_vector_[0][i] * Avogadro_Constant * lambda[i] * TOXCI[i];
                    nuclide_string << " ";
                    if (n_vector_[size][i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][i] * Avogadro_Constant * lambda[i] * TOXCI[i];
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    toxicity_node->InsertEndChild(nuclide);
                } else {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    for ( int j = 0; j < size_tot; ++j) {
                        if(n_vector_[j][i] < cutoff) {
                            nuclide_string << 0.0 <<" ";
                        } else {
                            nuclide_string << n_vector_[j][i] * Avogadro_Constant * lambda[i] * TOXCI[i] << " ";
                        }
                    }
                    if (n_vector_[size_tot][i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][i] * Avogadro_Constant * lambda[i] * TOXCI[i];
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    toxicity_node->InsertEndChild(nuclide);
                }
            }

        }

        if (if_print_fission_rate_ == true) {
            XMLElement* prod_node = doc.NewElement("NeuProdRate");
            prod_node->SetAttribute("unit","s-1");
            nuclides->InsertEndChild(prod_node);

            XMLElement *prod_tot = doc.NewElement("total");
            if (print_mode_ == 0) {
                stringstream prod_string;
                prod_string << scientific << setprecision(output_precision_);
                prod_string << prod_neu[0];
                prod_string << " ";
                prod_string << prod_neu.back();

                string prod_s(prod_string.str());

                XMLText *prod_text = doc.NewText(prod_s.c_str());
                prod_tot->InsertEndChild(prod_text);
                prod_node->InsertEndChild(prod_tot);
            } else {
                stringstream prod_string;
                prod_string << scientific << setprecision(output_precision_);
                for (int j = 0; j < size_tot; ++j) {
                    prod_string << prod_neu[j] << " ";
                }
                prod_string << prod_neu.back();

                string prod_s(prod_string.str());
                XMLText *prod_text = doc.NewText(prod_s.c_str());
                prod_tot->InsertEndChild(prod_text);
                prod_node->InsertEndChild(prod_tot);
            }

            for (int i = 0; i < nucl_size; ++ i) {
                XMLElement *nuclide = doc.NewElement("nuclide");
                int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                int NZ = nucl_id / 10000;
                int NA = (nucl_id - NZ * 10000) / 10;
                int NG = nucl_id - nucl_id / 10 * 10;
                string nucl_name;
                nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
                if(NG == 1) {
                    nucl_name += "m";
                }

                nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
                nuclide->SetAttribute("name",nucl_name.c_str());
                if(print_mode_ == 0) {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    nuclide_string << n_vector_[0][i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24 * flux_vector_[0];
                    nuclide_string << " ";
                    if (n_vector_[size][i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24 * flux_vector_[size_tot];
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    prod_node->InsertEndChild(nuclide);
                } else {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    for ( int j = 0; j < size_tot; ++j) {
                        if(n_vector_[j][i] < cutoff) {
                            nuclide_string << 0.0 <<" ";
                        } else {
                            nuclide_string << n_vector_[j][i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24 * flux_vector_[j] << " ";
                        }
                    }
                    if (n_vector_[size_tot][i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24 * flux_vector_[size_tot];
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    prod_node->InsertEndChild(nuclide);
                }
            }

        }


        if( if_print_absorption_rate_ == true) {
            XMLElement* abs_node = doc.NewElement("NeuAbsRate");
            abs_node->SetAttribute("unit","s-1");
            nuclides->InsertEndChild(abs_node);

            XMLElement *abs_tot = doc.NewElement("total");
            if(print_mode_ == 0) {
                stringstream abs_string;
                abs_string << scientific << setprecision(output_precision_);
                abs_string << abs_neu[0];
                abs_string << " ";
                abs_string << abs_neu.back();

                string abs_s(abs_string.str());

                XMLText * abs_text = doc.NewText(abs_s.c_str());
                abs_tot->InsertEndChild(abs_text);
                abs_node->InsertEndChild(abs_tot);
            } else {
                stringstream abs_string;
                abs_string << scientific << setprecision(output_precision_);
                for(int j = 0; j < size_tot; ++j) {
                    abs_string << abs_neu[j] << " ";
                }
                abs_string << abs_neu.back();

                string abs_s(abs_string.str());
                XMLText *abs_text = doc.NewText(abs_s.c_str());
                abs_tot->InsertEndChild(abs_text);
                abs_node->InsertEndChild(abs_tot);
            }

            for (int i = 0; i < nucl_size; ++ i) {
                XMLElement *nuclide = doc.NewElement("nuclide");
                int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                int NZ = nucl_id / 10000;
                int NA = (nucl_id - NZ * 10000) / 10;
                int NG = nucl_id - nucl_id / 10 * 10;
                string nucl_name;
                nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
                if(NG == 1) {
                    nucl_name += "m";
                }

                nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
                nuclide->SetAttribute("name",nucl_name.c_str());
                if(print_mode_ == 0) {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    nuclide_string << n_vector_[0][i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24 * flux_vector_[0];
                    nuclide_string << " ";
                    if (n_vector_[size][i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24 * flux_vector_[size_tot];
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    abs_node->InsertEndChild(nuclide);
                } else {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    for ( int j = 0; j < size_tot; ++j) {
                        if(n_vector_[j][i] < cutoff) {
                            nuclide_string << 0.0 <<" ";
                        } else {
                            nuclide_string << n_vector_[j][i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24 * flux_vector_[j] << " ";
                        }
                    }
                    if (n_vector_[size_tot][i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24 * flux_vector_[size_tot];
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    abs_node->InsertEndChild(nuclide);
                }
            }
        }
        if (if_tracking_stockage == true) {
            XMLElement *nuclides = doc.NewElement("Nuclides");
            nuclides->SetAttribute("zone","stockage");
            root->InsertEndChild(nuclides);

            XMLElement *concentration = doc.NewElement("Concentration");
            concentration->SetAttribute("unit",dens_unit_.c_str());
            nuclides->InsertEndChild(concentration);

            double convert_coeff;

            for ( int i = 0; i < nucl_size; ++i) {
                XMLElement *nuclide = doc.NewElement("nuclide");
                int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                int NZ = nucl_id / 10000;
                int NA = (nucl_id - NZ * 10000) / 10;
                int NG = nucl_id - nucl_id / 10 * 10;
                string nucl_name;
                nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
                if(NG == 1) {
                    nucl_name += "m";
                }

                if (dens_unit_ == "mol") {
                    convert_coeff = 1.0;
                } else if (dens_unit_ == "g") {
                    convert_coeff = double(NA);
                } else if (dens_unit_ == "kg") {
                    convert_coeff = double(NA) / 1000.0;
                } else if (dens_unit_ == "atom") {
                    convert_coeff = Avogadro_Constant;
                } else if (dens_unit_ == "atom/(barn-cm)") {
                    convert_coeff = Avogadro_Constant * 1.0e-24;
                }

                nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
                nuclide->SetAttribute("name",nucl_name.c_str());
                if (print_mode_ == 0) {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);
                    nuclide_string << n_vector_[0][i + nucl_size] * convert_coeff;
                    nuclide_string << " ";
                    if(n_vector_[size_tot][i + nucl_size] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][i + nucl_size] * convert_coeff;
                    }
                    string nuclide_s(nuclide_string.str());
                    XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    concentration->InsertEndChild(nuclide);
                } else {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);
                    for (int j = 0; j < size_tot; ++j) {
                        if (n_vector_[j][i + nucl_size] < cutoff) {
                            nuclide_string << 0.0 << " ";
                        } else {
                            nuclide_string << n_vector_[j][i + nucl_size] * convert_coeff << " ";
                        }
                    }
                    if (n_vector_[size_tot][i + nucl_size] < cutoff) {
                        nuclide_string << 0.0 ;
                    } else {

                        nuclide_string << n_vector_[size_tot][i + nucl_size] * convert_coeff;
                    }
                    string nuclide_s(nuclide_string.str());
                    XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    concentration->InsertEndChild(nuclide);
                }

            }
            if (if_print_activity_ == true) {
                XMLElement * act_node = doc.NewElement("RadioActivity");
                act_node->SetAttribute("unit","Ci");
                nuclides->InsertEndChild(act_node);

                vector< double > lambda;
                lambda = ModecNuclideLibrary.nuclide_library_vector_[1];

                for (int i = 0; i < nucl_size; ++ i) {
                    XMLElement *nuclide = doc.NewElement("nuclide");
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string nucl_name;
                    nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
                    if(NG == 1) {
                        nucl_name += "m";
                    }

                    nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
                    nuclide->SetAttribute("name",nucl_name.c_str());
                    if(print_mode_ == 0) {
                        stringstream nuclide_string;
                        nuclide_string << scientific << setprecision(output_precision_);

                        nuclide_string << n_vector_[0][i + nucl_size] * Avogadro_Constant * lambda[i] / (3.7e10);
                        nuclide_string << " ";
                        if (n_vector_[size][i + nucl_size] < cutoff) {
                            nuclide_string << 0.0;
                        } else {
                            nuclide_string << n_vector_[size_tot][i + nucl_size] * Avogadro_Constant * lambda[i] / (3.7e10);
                        }

                        string nuclide_s(nuclide_string.str());
                        XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                        nuclide->InsertEndChild(nuclide_text);
                        act_node->InsertEndChild(nuclide);
                    } else {
                        stringstream nuclide_string;
                        nuclide_string << scientific << setprecision(output_precision_);

                        for ( int j = 0; j < size_tot; ++j) {
                            if(n_vector_[j][i + nucl_size] < cutoff) {
                                nuclide_string << 0.0 <<" ";
                            } else {
                                nuclide_string << n_vector_[j][i + nucl_size] * Avogadro_Constant * lambda[i] / (3.7e10) << " ";
                            }
                        }
                        if (n_vector_[size_tot][i + nucl_size] < cutoff) {
                            nuclide_string << 0.0;
                        } else {
                            nuclide_string << n_vector_[size_tot][i + nucl_size] * Avogadro_Constant * lambda[i] / (3.7e10);
                        }

                        string nuclide_s(nuclide_string.str());
                        XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                        nuclide->InsertEndChild(nuclide_text);
                        act_node->InsertEndChild(nuclide);
                    }
                }

            }

            if (if_print_decayenergy_ == true) {

                XMLElement * decay_node = doc.NewElement("DecayHeat");
                decay_node->SetAttribute("unit","Watts");
                nuclides->InsertEndChild(decay_node);

                vector< double > lambda;
                lambda = ModecNuclideLibrary.nuclide_library_vector_[1];

                vector< double > Q;
                Q = ModecNuclideLibrary.nuclide_library_vector_[3];

                for (int i = 0; i < nucl_size; ++ i) {
                    XMLElement *nuclide = doc.NewElement("nuclide");
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string nucl_name;
                    nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
                    if(NG == 1) {
                        nucl_name += "m";
                    }

                    nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
                    nuclide->SetAttribute("name",nucl_name.c_str());
                    if(print_mode_ == 0) {
                        stringstream nuclide_string;
                        nuclide_string << scientific << setprecision(output_precision_);

                        nuclide_string << n_vector_[0][i + nucl_size] * Avogadro_Constant * lambda[i] * Q[i] * Electron_Coulomb * 1.0e6;
                        nuclide_string << " ";
                        if (n_vector_[size][i + nucl_size] < cutoff) {
                            nuclide_string << 0.0;
                        } else {
                            nuclide_string << n_vector_[size_tot][i + nucl_size] * Avogadro_Constant * lambda[i] * Q[i] * Electron_Coulomb * 1.0e6;
                        }

                        string nuclide_s(nuclide_string.str());
                        XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                        nuclide->InsertEndChild(nuclide_text);
                        decay_node->InsertEndChild(nuclide);
                    } else {
                        stringstream nuclide_string;
                        nuclide_string << scientific << setprecision(output_precision_);

                        for ( int j = 0; j < size_tot; ++j) {
                            if(n_vector_[j][i + nucl_size] < cutoff) {
                                nuclide_string << 0.0 <<" ";
                            } else {
                                nuclide_string << n_vector_[j][i + nucl_size] * Avogadro_Constant * lambda[i] * Q[i] * Electron_Coulomb * 1.0e6 << " ";
                            }
                        }
                        if (n_vector_[size_tot][i + nucl_size] < cutoff) {
                            nuclide_string << 0.0;
                        } else {
                            nuclide_string << n_vector_[size_tot][i + nucl_size] * Avogadro_Constant * lambda[i] * Q[i] * Electron_Coulomb * 1.0e6;
                        }

                        string nuclide_s(nuclide_string.str());
                        XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                        nuclide->InsertEndChild(nuclide_text);
                        decay_node->InsertEndChild(nuclide);
                    }
                }

            }

            if (if_print_ampc_ == true) {

                XMLElement * ampc_node = doc.NewElement("AMPC");
                ampc_node->SetAttribute("unit","m3-air/Bq");
                nuclides->InsertEndChild(ampc_node);

                vector< double > lambda;
                lambda = ModecNuclideLibrary.nuclide_library_vector_[1];

                vector< double > AMPC;
                AMPC = ModecNuclideLibrary.nuclide_library_vector_[4];

                for (int i = 0; i < nucl_size; ++ i) {
                    XMLElement *nuclide = doc.NewElement("nuclide");
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string nucl_name;
                    nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
                    if(NG == 1) {
                        nucl_name += "m";
                    }

                    nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
                    nuclide->SetAttribute("name",nucl_name.c_str());
                    if(print_mode_ == 0) {
                        stringstream nuclide_string;
                        nuclide_string << scientific << setprecision(output_precision_);
                        nuclide_string << n_vector_[0][i + nucl_size] * Avogadro_Constant * lambda[i] / 3.7e10 / AMPC[i];
                        nuclide_string << " ";
                        if (n_vector_[size][i + nucl_size] < cutoff) {
                            nuclide_string << 0.0;
                        } else {
                            nuclide_string << n_vector_[size_tot][i + nucl_size] * Avogadro_Constant * lambda[i] / 3.7e10 / AMPC[i];
                        }

                        string nuclide_s(nuclide_string.str());
                        XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                        nuclide->InsertEndChild(nuclide_text);
                        ampc_node->InsertEndChild(nuclide);
                    } else {
                        stringstream nuclide_string;
                        nuclide_string << scientific << setprecision(output_precision_);

                        for ( int j = 0; j < size_tot; ++j) {
                            if(n_vector_[j][i + nucl_size] < cutoff) {
                                nuclide_string << 0.0 <<" ";
                            } else {
                                nuclide_string << n_vector_[j][i + nucl_size] * Avogadro_Constant * lambda[i] / 3.7e10 / AMPC[i] << " ";
                            }
                        }
                        if (n_vector_[size_tot][i + nucl_size] < cutoff) {
                            nuclide_string << 0.0;
                        } else {
                            nuclide_string << n_vector_[size_tot][i + nucl_size] * Avogadro_Constant * lambda[i] / 3.7e10 / AMPC[i];
                        }

                        string nuclide_s(nuclide_string.str());
                        XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                        nuclide->InsertEndChild(nuclide_text);
                        ampc_node->InsertEndChild(nuclide);
                    }
                }

            }

            if (if_print_wmpc_ == true) {

                XMLElement * wmpc_node = doc.NewElement("WMPC");
                wmpc_node->SetAttribute("unit","m3-water/Bq");
                nuclides->InsertEndChild(wmpc_node);

                vector< double > lambda;
                lambda = ModecNuclideLibrary.nuclide_library_vector_[1];

                vector< double > WMPC;
                WMPC = ModecNuclideLibrary.nuclide_library_vector_[5];

                for (int i = 0; i < nucl_size; ++ i) {
                    XMLElement *nuclide = doc.NewElement("nuclide");
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string nucl_name;
                    nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
                    if(NG == 1) {
                        nucl_name += "m";
                    }

                    nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
                    nuclide->SetAttribute("name",nucl_name.c_str());
                    if(print_mode_ == 0) {
                        stringstream nuclide_string;
                        nuclide_string << scientific << setprecision(output_precision_);

                        nuclide_string << n_vector_[0][i] * Avogadro_Constant * lambda[i] / 3.7e10 / WMPC[i];
                        nuclide_string << " ";
                        if (n_vector_[size][i + nucl_size] < cutoff) {
                            nuclide_string << 0.0;
                        } else {
                            nuclide_string << n_vector_[size_tot][i + nucl_size] * Avogadro_Constant * lambda[i] / 3.7e10 / WMPC[i];
                        }

                        string nuclide_s(nuclide_string.str());
                        XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                        nuclide->InsertEndChild(nuclide_text);
                        wmpc_node->InsertEndChild(nuclide);
                    } else {
                        stringstream nuclide_string;
                        nuclide_string << scientific << setprecision(output_precision_);

                        for ( int j = 0; j < size_tot; ++j) {
                            if(n_vector_[j][i + nucl_size] < cutoff) {
                                nuclide_string << 0.0 <<" ";
                            } else {
                                nuclide_string << n_vector_[j][i + nucl_size] * Avogadro_Constant * lambda[i] / 3.7e10 / WMPC[i] << " ";
                            }
                        }
                        if (n_vector_[size_tot][i + nucl_size] < cutoff) {
                            nuclide_string << 0.0;
                        } else {
                            nuclide_string << n_vector_[size_tot][i + nucl_size] * Avogadro_Constant * lambda[i] / 3.7e10 / WMPC[i];
                        }

                        string nuclide_s(nuclide_string.str());
                        XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                        nuclide->InsertEndChild(nuclide_text);
                        wmpc_node->InsertEndChild(nuclide);
                    }
                }

            }

            if (if_print_toxicity_ == true) {

                XMLElement * toxicity_node = doc.NewElement("RadioToxicity");
                toxicity_node->SetAttribute("unit","Sv");
                nuclides->InsertEndChild(toxicity_node);

                vector< double > lambda;
                lambda = ModecNuclideLibrary.nuclide_library_vector_[1];

                vector< double > TOXCI;
                TOXCI = ModecNuclideLibrary.nuclide_library_vector_[10];

                for (int i = 0; i < nucl_size; ++ i) {
                    XMLElement *nuclide = doc.NewElement("nuclide");
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string nucl_name;
                    nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
                    if(NG == 1) {
                        nucl_name += "m";
                    }

                    nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
                    nuclide->SetAttribute("name",nucl_name.c_str());
                    if(print_mode_ == 0) {
                        stringstream nuclide_string;
                        nuclide_string << scientific << setprecision(output_precision_);

                        nuclide_string << n_vector_[0][i + nucl_size] * Avogadro_Constant * lambda[i] * TOXCI[i];
                        nuclide_string << " ";
                        if (n_vector_[size][i + nucl_size] < cutoff) {
                            nuclide_string << 0.0;
                        } else {
                            nuclide_string << n_vector_[size_tot][i + nucl_size] * Avogadro_Constant * lambda[i] * TOXCI[i];
                        }

                        string nuclide_s(nuclide_string.str());
                        XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                        nuclide->InsertEndChild(nuclide_text);
                        toxicity_node->InsertEndChild(nuclide);
                    } else {
                        stringstream nuclide_string;
                        nuclide_string << scientific << setprecision(output_precision_);

                        for ( int j = 0; j < size_tot; ++j) {
                            if(n_vector_[j][i + nucl_size] < cutoff) {
                                nuclide_string << 0.0 <<" ";
                            } else {
                                nuclide_string << n_vector_[j][i + nucl_size] * Avogadro_Constant * lambda[i] * TOXCI[i] << " ";
                            }
                        }
                        if (n_vector_[size_tot][i + nucl_size] < cutoff) {
                            nuclide_string << 0.0;
                        } else {
                            nuclide_string << n_vector_[size_tot][i + nucl_size] * Avogadro_Constant * lambda[i] * TOXCI[i];
                        }

                        string nuclide_s(nuclide_string.str());
                        XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                        nuclide->InsertEndChild(nuclide_text);
                        toxicity_node->InsertEndChild(nuclide);
                    }
                }

            }

        }
    } else if (if_flow_mode_ == 1) {
        XMLElement *time = doc.NewElement("Time");
        time->SetAttribute("unit","s");
        if(print_mode_ == 0) {
            stringstream time_string;
            time_string << scientific << setprecision(output_precision_);
            time_string << cum_time[0];
            time_string << " ";
            time_string << cum_time.back();

            string time_s(time_string.str());
            XMLText *time_text = doc.NewText(time_s.c_str());
            time->InsertEndChild(time_text);
            root->InsertEndChild(time);
        } else {
            stringstream time_string;
            time_string << scientific << setprecision(output_precision_);
            for (int j = 0; j < size_tot; ++j) {
                time_string << cum_time[j] << " ";
            }
            time_string << cum_time.back();

            string time_s(time_string.str());
            //time_string >> time_s;
            XMLText *time_text = doc.NewText(time_s.c_str());
            time->InsertEndChild(time_text);
            root->InsertEndChild(time);
        }

        XMLElement *power = doc.NewElement("Power");
        power->SetAttribute("unit","MW");
        if(print_mode_ == 0) {
            stringstream power_string;
            power_string << scientific << setprecision(output_precision_);
            power_string << power_vector_[0];
            power_string << " ";
            power_string << power_vector_.back();

            string power_s(power_string.str());
            XMLText * power_text = doc.NewText(power_s.c_str());
            power->InsertEndChild(power_text);
            root->InsertEndChild(power);
        } else {
            stringstream power_string;
            power_string << scientific << setprecision(output_precision_);
            for (int j = 0; j < size_tot; ++j) {
                power_string << power_vector_[j] << " ";
            }
            power_string << power_vector_.back();

            string power_s(power_string.str());
            XMLText *power_text = doc.NewText(power_s.c_str());
            power->InsertEndChild(power_text);
            root->InsertEndChild(power);
        }

        XMLElement *flux = doc.NewElement("Flux");
        flux->SetAttribute("unit","cm-2*s-1");
        if(print_mode_ == 0) {
            stringstream flux_string;
            flux_string << scientific << setprecision(output_precision_);
            flux_string << flux_vector_[0];
            flux_string << " ";
            flux_string << flux_vector_.back();

            string flux_s(flux_string.str());
            XMLText * flux_text = doc.NewText(flux_s.c_str());
            flux->InsertEndChild(flux_text);
            root->InsertEndChild(flux);
        } else {
            stringstream flux_string;
            flux_string << scientific << setprecision(output_precision_);
            for (int j = 0; j < size_tot; ++j) {
                flux_string << flux_vector_[j] << " ";
            }
            flux_string << flux_vector_.back();

            string flux_s(flux_string.str());
            XMLText *flux_text = doc.NewText(flux_s.c_str());
            flux->InsertEndChild(flux_text);
            root->InsertEndChild(flux);
        }

        vector<double> kinf;
        vector<double> prod_neu;
        vector<double> abs_neu;
        if (if_print_kinf_ != 0 || if_print_fission_rate_ != 0 || if_print_absorption_rate_ != 0) {
            kinf.resize(size_tot + 1);
            prod_neu.resize(size_tot + 1);
            abs_neu.resize(size_tot + 1);

            for (int j = 0; j <= size_tot; ++j) {
                for (int i = 0; i < nucl_size; ++i) {
                    prod_neu[j] += n_vector_[j][i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24;
                    abs_neu[j] += n_vector_[j][i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24;
                }

                prod_neu[j] *= flux_vector_[j];
                abs_neu[j] *= flux_vector_[j];

                if (abs_neu[j] != 0.0) {
                    kinf[j] = prod_neu[j] / abs_neu[j];
                }
            }

            if (if_print_kinf_ == true) {
                XMLElement* kinf_node = doc.NewElement("K-infinite");
                if (print_mode_ == 0) {
                    stringstream kinf_string;
                    kinf_string << scientific << setprecision(output_precision_);
                    kinf_string << kinf[0];
                    kinf_string << " ";
                    kinf_string << kinf.back();

                    string kinf_s(kinf_string.str());

                    XMLText *kinf_text = doc.NewText(kinf_s.c_str());
                    kinf_node->InsertEndChild(kinf_text);
                    root->InsertEndChild(kinf_node);
                } else {
                    stringstream kinf_string;
                    kinf_string << scientific << setprecision(output_precision_);
                    for (int j = 0; j < size_tot; ++j) {
                        kinf_string << kinf[j] << " ";
                    }
                    kinf_string << kinf.back();

                    string kinf_s(kinf_string.str());
                    XMLText *kinf_text = doc.NewText(kinf_s.c_str());
                    kinf_node->InsertEndChild(kinf_text);
                    root->InsertEndChild(kinf_node);
                }
            }

        }

        XMLElement *nuclides = doc.NewElement("Nuclides");
        nuclides->SetAttribute("zone","core");
        root->InsertEndChild(nuclides);

        XMLElement *nuclides_loop = doc.NewElement("Nuclides");
        nuclides_loop->SetAttribute("zone","loop");
        root->InsertEndChild(nuclides_loop);

        XMLElement *concentration = doc.NewElement("Concentration");
        concentration->SetAttribute("unit",dens_unit_.c_str());
        nuclides->InsertEndChild(concentration);

        XMLElement *concentration_loop = doc.NewElement("Concentration");
        concentration_loop->SetAttribute("unit",dens_unit_.c_str());
        nuclides_loop->InsertEndChild(concentration_loop);

        double convert_coeff;

        for ( int i = 0; i < nucl_size; ++i) {
            XMLElement *nuclide = doc.NewElement("nuclide");
            XMLElement *nuclide_loop = doc.NewElement("nuclide");

            int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
            int NZ = nucl_id / 10000;
            int NA = (nucl_id - NZ * 10000) / 10;
            int NG = nucl_id - nucl_id / 10 * 10;
            string nucl_name;
            nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
            if(NG == 1) {
                nucl_name += "m";
            }

            if (dens_unit_ == "mol") {
                convert_coeff = 1.0;
            } else if (dens_unit_ == "g") {
                convert_coeff = double(NA);
            } else if (dens_unit_ == "kg") {
                convert_coeff = double(NA) / 1000.0;
            } else if (dens_unit_ == "atom") {
                convert_coeff = Avogadro_Constant;
            } else if (dens_unit_ == "atom/(barn-cm)") {
                convert_coeff = Avogadro_Constant * 1.0e-24;
            }

            nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
            nuclide->SetAttribute("name",nucl_name.c_str());

            nuclide_loop->SetAttribute("zai",to_string(nucl_id).c_str());
            nuclide_loop->SetAttribute("name",nucl_name.c_str());

            if (print_mode_ == 0) {
                stringstream nuclide_string;
                stringstream nuclide_loop_string;

                nuclide_string << scientific << setprecision(output_precision_);
                nuclide_string << n_vector_[0][2 * i] * convert_coeff;
                nuclide_string << " ";
                if(n_vector_[size_tot][2 * i] < cutoff) {
                    nuclide_string << 0.0;
                } else {
                    nuclide_string << n_vector_[size_tot][2 * i] * convert_coeff;
                }
                string nuclide_s(nuclide_string.str());
                XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                nuclide->InsertEndChild(nuclide_text);
                concentration->InsertEndChild(nuclide);

                nuclide_loop_string << scientific << setprecision(output_precision_);
                nuclide_loop_string << n_vector_[0][2 * i + 1] * convert_coeff;
                nuclide_loop_string << " ";
                if(n_vector_[size_tot][2 * i + 1] < cutoff) {
                    nuclide_loop_string << 0.0;
                } else {
                    nuclide_loop_string << n_vector_[size_tot][2 * i + 1] * convert_coeff;
                }
                string nuclide_loop_s(nuclide_loop_string.str());
                XMLText *nuclide_loop_text = doc.NewText(nuclide_loop_s.c_str());
                nuclide_loop->InsertEndChild(nuclide_loop_text);
                concentration_loop->InsertEndChild(nuclide_loop);
            } else {
                stringstream nuclide_string;
                stringstream nuclide_loop_string;

                nuclide_string << scientific << setprecision(output_precision_);
                for (int j = 0; j < size_tot; ++j) {
                    if (n_vector_[j][2 * i] < cutoff) {
                        nuclide_string << 0.0 << " ";
                    } else {
                        nuclide_string << n_vector_[j][2 * i] * convert_coeff << " ";
                    }
                }
                if (n_vector_[size_tot][2 * i] < cutoff) {
                    nuclide_string << 0.0 ;
                } else {

                    nuclide_string << n_vector_[size_tot][2 * i] * convert_coeff;
                }
                string nuclide_s(nuclide_string.str());
                XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                nuclide->InsertEndChild(nuclide_text);
                concentration->InsertEndChild(nuclide);

                /// 外回路核素浓度
                nuclide_loop_string << scientific << setprecision(output_precision_);
                for (int j = 0; j < size_tot; ++j) {
                    if (n_vector_[j][2 * i + 1] < cutoff) {
                        nuclide_loop_string << 0.0 << " ";
                    } else {
                        nuclide_loop_string << n_vector_[j][2 * i + 1] * convert_coeff << " ";
                    }
                }
                if (n_vector_[size_tot][2 * i + 1] < cutoff) {
                    nuclide_loop_string << 0.0 ;
                } else {

                    nuclide_loop_string << n_vector_[size_tot][2 * i + 1] * convert_coeff;
                }
                string nuclide_loop_s(nuclide_loop_string.str());
                XMLText * nuclide_loop_text = doc.NewText(nuclide_loop_s.c_str());
                nuclide_loop->InsertEndChild(nuclide_loop_text);
                concentration_loop->InsertEndChild(nuclide_loop);
            }

        }
        if (if_print_activity_ == true) {
            XMLElement * act_node = doc.NewElement("RadioActivity");
            XMLElement * act_loop_node = doc.NewElement("RadioActivity");

            act_node->SetAttribute("unit","Ci");
            nuclides->InsertEndChild(act_node);

            /// 外回路
            act_loop_node->SetAttribute("unit","Ci");
            nuclides_loop->InsertEndChild(act_loop_node);

            vector< double > lambda;
            lambda = ModecNuclideLibrary.nuclide_library_vector_[1];

            for (int i = 0; i < nucl_size; ++ i) {
                XMLElement *nuclide = doc.NewElement("nuclide");
                XMLElement *nuclide_loop = doc.NewElement("nuclide");

                int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                int NZ = nucl_id / 10000;
                int NA = (nucl_id - NZ * 10000) / 10;
                int NG = nucl_id - nucl_id / 10 * 10;
                string nucl_name;
                nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
                if(NG == 1) {
                    nucl_name += "m";
                }

                nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
                nuclide->SetAttribute("name",nucl_name.c_str());

                nuclide_loop->SetAttribute("zai",to_string(nucl_id).c_str());
                nuclide_loop->SetAttribute("name",nucl_name.c_str());

                if(print_mode_ == 0) {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    nuclide_string << n_vector_[0][2 * i] * Avogadro_Constant * lambda[i] / (3.7e10);
                    nuclide_string << " ";
                    if (n_vector_[size][2 * i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][2 * i] * Avogadro_Constant * lambda[i] / (3.7e10);
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    act_node->InsertEndChild(nuclide);

                    /// 外回路
                    stringstream nuclide_loop_string;
                    nuclide_loop_string << scientific << setprecision(output_precision_);

                    nuclide_loop_string << n_vector_[0][2 * i + 1] * Avogadro_Constant * lambda[i] / (3.7e10);
                    nuclide_loop_string << " ";
                    if (n_vector_[size][2 * i + 1] < cutoff) {
                        nuclide_loop_string << 0.0;
                    } else {
                        nuclide_loop_string << n_vector_[size_tot][2 * i + 1] * Avogadro_Constant * lambda[i] / (3.7e10);
                    }

                    string nuclide_loop_s(nuclide_loop_string.str());
                    XMLText * nuclide_loop_text = doc.NewText(nuclide_loop_s.c_str());
                    nuclide_loop->InsertEndChild(nuclide_loop_text);
                    act_loop_node->InsertEndChild(nuclide_loop);
                } else {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    for ( int j = 0; j < size_tot; ++j) {
                        if(n_vector_[j][2 * i] < cutoff) {
                            nuclide_string << 0.0 <<" ";
                        } else {
                            nuclide_string << n_vector_[j][2 * i] * Avogadro_Constant * lambda[i] / (3.7e10) << " ";
                        }
                    }
                    if(n_vector_[size_tot][2 * i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][2 * i] * Avogadro_Constant * lambda[i] / (3.7e10);
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    act_node->InsertEndChild(nuclide);

                    /// 外回路
                    stringstream nuclide_loop_string;
                    nuclide_loop_string << scientific << setprecision(output_precision_);

                    for ( int j = 0; j < size_tot; ++j) {
                        if(n_vector_[j][2 * i + 1] < cutoff) {
                            nuclide_loop_string << 0.0 <<" ";
                        } else {
                            nuclide_loop_string << n_vector_[j][2 * i + 1] * Avogadro_Constant * lambda[i] / (3.7e10) << " ";
                        }
                    }
                    if(n_vector_[size_tot][2 * i + 1] < cutoff) {
                        nuclide_loop_string << 0.0;
                    } else {
                        nuclide_loop_string << n_vector_[size_tot][2 * i + 1] * Avogadro_Constant * lambda[i] / (3.7e10);
                    }

                    string nuclide_loop_s(nuclide_loop_string.str());
                    XMLText *nuclide_loop_text = doc.NewText(nuclide_loop_s.c_str());
                    nuclide_loop->InsertEndChild(nuclide_loop_text);
                    act_loop_node->InsertEndChild(nuclide_loop);
                }
            }

        }

        if (if_print_decayenergy_ == true) {

            XMLElement * decay_node = doc.NewElement("DecayHeat");
            decay_node->SetAttribute("unit","Watts");
            nuclides->InsertEndChild(decay_node);

            XMLElement * decay_loop_node = doc.NewElement("DecayHeat");
            decay_loop_node->SetAttribute("unit","Watts");
            nuclides_loop->InsertEndChild(decay_loop_node);

            vector< double > lambda;
            lambda = ModecNuclideLibrary.nuclide_library_vector_[1];

            vector< double > Q;
            Q = ModecNuclideLibrary.nuclide_library_vector_[3];

            for (int i = 0; i < nucl_size; ++ i) {
                XMLElement *nuclide = doc.NewElement("nuclide");
                XMLElement *nuclide_loop = doc.NewElement("nuclide");

                int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                int NZ = nucl_id / 10000;
                int NA = (nucl_id - NZ * 10000) / 10;
                int NG = nucl_id - nucl_id / 10 * 10;
                string nucl_name;
                nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
                if(NG == 1) {
                    nucl_name += "m";
                }

                nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
                nuclide->SetAttribute("name",nucl_name.c_str());

                nuclide_loop->SetAttribute("zai",to_string(nucl_id).c_str());
                nuclide_loop->SetAttribute("name",nucl_name.c_str());

                if(print_mode_ == 0) {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    nuclide_string << n_vector_[0][2 * i] * Avogadro_Constant * lambda[i] * Q[i] * Electron_Coulomb * 1.0e6;
                    nuclide_string << " ";
                    if (n_vector_[size][2 * i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][2 * i] * Avogadro_Constant * lambda[i] * Q[i] * Electron_Coulomb * 1.0e6;
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    decay_node->InsertEndChild(nuclide);

                    /// 外回路
                    stringstream nuclide_loop_string;
                    nuclide_loop_string << scientific << setprecision(output_precision_);

                    nuclide_loop_string << n_vector_[0][2 * i + 1] * Avogadro_Constant * lambda[i] * Q[i] * Electron_Coulomb * 1.0e6;
                    nuclide_loop_string << " ";
                    if (n_vector_[size][2 * i + 1] < cutoff) {
                        nuclide_loop_string << 0.0;
                    } else {
                        nuclide_loop_string << n_vector_[size_tot][2 * i + 1] * Avogadro_Constant * lambda[i] * Q[i] * Electron_Coulomb * 1.0e6;
                    }

                    string nuclide_loop_s(nuclide_loop_string.str());
                    XMLText * nuclide_loop_text = doc.NewText(nuclide_loop_s.c_str());
                    nuclide_loop->InsertEndChild(nuclide_loop_text);
                    decay_loop_node->InsertEndChild(nuclide_loop);
                } else {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    for ( int j = 0; j < size_tot; ++j) {
                        if(n_vector_[j][2 * i] < cutoff) {
                            nuclide_string << 0.0 <<" ";
                        } else {
                            nuclide_string << n_vector_[j][2 * i] * Avogadro_Constant * lambda[i] * Q[i] * Electron_Coulomb * 1.0e6 << " ";
                        }
                    }
                    if(n_vector_[size_tot][2 * i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][2 * i] * Avogadro_Constant * lambda[i] * Q[i] * Electron_Coulomb * 1.0e6;
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    decay_node->InsertEndChild(nuclide);

                    /// 外回路
                    stringstream nuclide_loop_string;
                    nuclide_loop_string << scientific << setprecision(output_precision_);

                    for ( int j = 0; j < size_tot; ++j) {
                        if(n_vector_[j][2 * i + 1] < cutoff) {
                            nuclide_loop_string << 0.0 <<" ";
                        } else {
                            nuclide_loop_string << n_vector_[j][2 * i + 1] * Avogadro_Constant * lambda[i] * Q[i] * Electron_Coulomb * 1.0e6 << " ";
                        }
                    }
                    if( n_vector_[size_tot][2 * i + 1] < cutoff) {
                        nuclide_loop_string << 0.0;
                    } else {
                        nuclide_loop_string << n_vector_[size_tot][2 * i + 1] * Avogadro_Constant * lambda[i] * Q[i] * Electron_Coulomb * 1.0e6;
                    }

                    string nuclide_loop_s(nuclide_loop_string.str());
                    XMLText *nuclide_loop_text = doc.NewText(nuclide_loop_s.c_str());
                    nuclide_loop->InsertEndChild(nuclide_loop_text);
                    decay_loop_node->InsertEndChild(nuclide_loop);
                }
            }

        }

        if (if_print_ampc_ == true) {

            XMLElement * ampc_node = doc.NewElement("AMPC");
            ampc_node->SetAttribute("unit","m3-air/Bq");
            nuclides->InsertEndChild(ampc_node);

            XMLElement * ampc_loop_node = doc.NewElement("AMPC");
            ampc_loop_node->SetAttribute("unit","m3-air/Bq");
            nuclides_loop->InsertEndChild(ampc_loop_node);

            vector< double > lambda;
            lambda = ModecNuclideLibrary.nuclide_library_vector_[1];

            vector< double > AMPC;
            AMPC = ModecNuclideLibrary.nuclide_library_vector_[4];

            for (int i = 0; i < nucl_size; ++ i) {
                XMLElement *nuclide = doc.NewElement("nuclide");
                XMLElement *nuclide_loop = doc.NewElement("nuclide");

                int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                int NZ = nucl_id / 10000;
                int NA = (nucl_id - NZ * 10000) / 10;
                int NG = nucl_id - nucl_id / 10 * 10;
                string nucl_name;
                nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
                if(NG == 1) {
                    nucl_name += "m";
                }

                nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
                nuclide->SetAttribute("name",nucl_name.c_str());

                nuclide_loop->SetAttribute("zai",to_string(nucl_id).c_str());
                nuclide_loop->SetAttribute("name",nucl_name.c_str());

                if(print_mode_ == 0) {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);
                    nuclide_string << n_vector_[0][2 * i] * Avogadro_Constant * lambda[i] / 3.7e10 / AMPC[i];
                    nuclide_string << " ";
                    if (n_vector_[size][2 * i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][2 * i] * Avogadro_Constant * lambda[i] / 3.7e10 / AMPC[i];
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    ampc_node->InsertEndChild(nuclide);

                    /// 外回路
                    stringstream nuclide_loop_string;
                    nuclide_loop_string << scientific << setprecision(output_precision_);
                    nuclide_loop_string << n_vector_[0][2 * i + 1] * Avogadro_Constant * lambda[i] / 3.7e10 / AMPC[i];
                    nuclide_loop_string << " ";
                    if (n_vector_[size][2 * i + 1] < cutoff) {
                        nuclide_loop_string << 0.0;
                    } else {
                        nuclide_loop_string << n_vector_[size_tot][2 * i + 1] * Avogadro_Constant * lambda[i] / 3.7e10 / AMPC[i];
                    }

                    string nuclide_loop_s(nuclide_loop_string.str());
                    XMLText * nuclide_loop_text = doc.NewText(nuclide_loop_s.c_str());
                    nuclide_loop->InsertEndChild(nuclide_loop_text);
                    ampc_loop_node->InsertEndChild(nuclide_loop);
                } else {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    for ( int j = 0; j < size_tot; ++j) {
                        if(n_vector_[j][2 * i] < cutoff) {
                            nuclide_string << 0.0 <<" ";
                        } else {
                            nuclide_string << n_vector_[j][2 * i] * Avogadro_Constant * lambda[i] / 3.7e10 / AMPC[i] << " ";
                        }
                    }
                    if (n_vector_[size_tot][2 * i + 1] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][2 * i] * Avogadro_Constant * lambda[i] / 3.7e10 / AMPC[i];
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    ampc_node->InsertEndChild(nuclide);

                    /// 外回路
                    stringstream nuclide_loop_string;
                    nuclide_loop_string << scientific << setprecision(output_precision_);

                    for ( int j = 0; j < size_tot; ++j) {
                        if(n_vector_[j][2 * i + 1] < cutoff) {
                            nuclide_loop_string << 0.0 <<" ";
                        } else {
                            nuclide_loop_string << n_vector_[j][2 * i + 1] * Avogadro_Constant * lambda[i] / 3.7e10 / AMPC[i] << " ";
                        }
                    }
                    if (n_vector_[size_tot][2 * i + 1] < cutoff) {
                        nuclide_loop_string << 0.0;
                    } else {
                        nuclide_loop_string << n_vector_[size_tot][2 * i + 1] * Avogadro_Constant * lambda[i] / 3.7e10 / AMPC[i];
                    }

                    string nuclide_loop_s(nuclide_loop_string.str());
                    XMLText *nuclide_loop_text = doc.NewText(nuclide_loop_s.c_str());
                    nuclide_loop->InsertEndChild(nuclide_loop_text);
                    ampc_loop_node->InsertEndChild(nuclide_loop);
                }
            }

        }

        if (if_print_wmpc_ == true) {

            XMLElement * wmpc_node = doc.NewElement("WMPC");
            XMLElement * wmpc_loop_node = doc.NewElement("WMPC");

            wmpc_node->SetAttribute("unit","m3-water/Bq");
            nuclides->InsertEndChild(wmpc_node);

            wmpc_loop_node->SetAttribute("unit","m3-water/Bq");
            nuclides_loop->InsertEndChild(wmpc_loop_node);

            vector< double > lambda;
            lambda = ModecNuclideLibrary.nuclide_library_vector_[1];

            vector< double > WMPC;
            WMPC = ModecNuclideLibrary.nuclide_library_vector_[5];

            for (int i = 0; i < nucl_size; ++ i) {
                XMLElement *nuclide = doc.NewElement("nuclide");
                XMLElement *nuclide_loop = doc.NewElement("nuclide");

                int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                int NZ = nucl_id / 10000;
                int NA = (nucl_id - NZ * 10000) / 10;
                int NG = nucl_id - nucl_id / 10 * 10;
                string nucl_name;
                nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
                if(NG == 1) {
                    nucl_name += "m";
                }

                nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
                nuclide->SetAttribute("name",nucl_name.c_str());

                nuclide_loop->SetAttribute("zai",to_string(nucl_id).c_str());
                nuclide_loop->SetAttribute("name",nucl_name.c_str());

                if(print_mode_ == 0) {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    nuclide_string << n_vector_[0][2 * i] * Avogadro_Constant * lambda[i] / 3.7e10 / WMPC[i];
                    nuclide_string << " ";
                    if (n_vector_[size][2 * i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][2 * i] * Avogadro_Constant * lambda[i] / 3.7e10 / WMPC[i];
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    wmpc_node->InsertEndChild(nuclide);

                    /// 外回路
                    stringstream nuclide_loop_string;
                    nuclide_loop_string << scientific << setprecision(output_precision_);

                    nuclide_loop_string << n_vector_[0][2 * i + 1] * Avogadro_Constant * lambda[i] / 3.7e10 / WMPC[i];
                    nuclide_loop_string << " ";
                    if (n_vector_[size][2 * i + 1] < cutoff) {
                        nuclide_loop_string << 0.0;
                    } else {
                        nuclide_loop_string << n_vector_[size_tot][2 * i + 1] * Avogadro_Constant * lambda[i] / 3.7e10 / WMPC[i];
                    }

                    string nuclide_loop_s(nuclide_loop_string.str());
                    XMLText * nuclide_loop_text = doc.NewText(nuclide_loop_s.c_str());
                    nuclide_loop->InsertEndChild(nuclide_loop_text);
                    wmpc_loop_node->InsertEndChild(nuclide_loop);
                } else {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    for ( int j = 0; j < size_tot; ++j) {
                        if(n_vector_[j][2 * i] < cutoff) {
                            nuclide_string << 0.0 <<" ";
                        } else {
                            nuclide_string << n_vector_[j][2 * i] * Avogadro_Constant * lambda[i] / 3.7e10 / WMPC[i] << " ";
                        }
                    }
                    if ( n_vector_[size_tot][2 * i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][2 * i] * Avogadro_Constant * lambda[i] / 3.7e10 / WMPC[i];
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    wmpc_node->InsertEndChild(nuclide);

                    /// 外回路
                    stringstream nuclide_loop_string;
                    nuclide_loop_string << scientific << setprecision(output_precision_);

                    for ( int j = 0; j < size_tot; ++j) {
                        if(n_vector_[j][2 * i + 1] < cutoff) {
                            nuclide_loop_string << 0.0 <<" ";
                        } else {
                            nuclide_loop_string << n_vector_[j][2 * i + 1] * Avogadro_Constant * lambda[i] / 3.7e10 / WMPC[i] << " ";
                        }
                    }
                    if ( n_vector_[size_tot][2 * i + 1] < cutoff) {
                        nuclide_loop_string << 0.0;
                    } else {
                        nuclide_loop_string << n_vector_[size_tot][2 * i + 1] * Avogadro_Constant * lambda[i] / 3.7e10 / WMPC[i];
                    }

                    string nuclide_loop_s(nuclide_loop_string.str());
                    XMLText *nuclide_loop_text = doc.NewText(nuclide_loop_s.c_str());
                    nuclide_loop->InsertEndChild(nuclide_loop_text);
                    wmpc_loop_node->InsertEndChild(nuclide_loop);
                }
            }

        }

        if (if_print_toxicity_ == true) {

            XMLElement * toxicity_node = doc.NewElement("RadioToxicity");
            XMLElement * toxicity_loop_node = doc.NewElement("RadioToxicity");

            toxicity_node->SetAttribute("unit","Sv");
            nuclides->InsertEndChild(toxicity_node);

            toxicity_loop_node->SetAttribute("unit","Sv");
            nuclides_loop->InsertEndChild(toxicity_loop_node);

            vector< double > lambda;
            lambda = ModecNuclideLibrary.nuclide_library_vector_[1];

            vector< double > TOXCI;
            TOXCI = ModecNuclideLibrary.nuclide_library_vector_[10];

            for (int i = 0; i < nucl_size; ++ i) {
                XMLElement *nuclide = doc.NewElement("nuclide");
                XMLElement *nuclide_loop = doc.NewElement("nuclide");

                int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                int NZ = nucl_id / 10000;
                int NA = (nucl_id - NZ * 10000) / 10;
                int NG = nucl_id - nucl_id / 10 * 10;
                string nucl_name;
                nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
                if(NG == 1) {
                    nucl_name += "m";
                }

                nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
                nuclide->SetAttribute("name",nucl_name.c_str());

                nuclide_loop->SetAttribute("zai",to_string(nucl_id).c_str());
                nuclide_loop->SetAttribute("name",nucl_name.c_str());

                if(print_mode_ == 0) {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    nuclide_string << n_vector_[0][2 * i] * Avogadro_Constant * lambda[i] * TOXCI[i];
                    nuclide_string << " ";
                    if (n_vector_[size][2 * i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][2 * i] * Avogadro_Constant * lambda[i] * TOXCI[i];
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    toxicity_node->InsertEndChild(nuclide);

                    /// 外回路
                    stringstream nuclide_loop_string;
                    nuclide_loop_string << scientific << setprecision(output_precision_);

                    nuclide_loop_string << n_vector_[0][2 * i + 1] * Avogadro_Constant * lambda[i] * TOXCI[i];
                    nuclide_loop_string << " ";
                    if (n_vector_[size][2 * i + 1] < cutoff) {
                        nuclide_loop_string << 0.0;
                    } else {
                        nuclide_loop_string << n_vector_[size_tot][2 * i + 1] * Avogadro_Constant * lambda[i] * TOXCI[i];
                    }

                    string nuclide_loop_s(nuclide_loop_string.str());
                    XMLText * nuclide_loop_text = doc.NewText(nuclide_loop_s.c_str());
                    nuclide_loop->InsertEndChild(nuclide_loop_text);
                    toxicity_loop_node->InsertEndChild(nuclide_loop);
                } else {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    for ( int j = 0; j < size_tot; ++j) {
                        if(n_vector_[j][2 * i] < cutoff) {
                            nuclide_string << 0.0 <<" ";
                        } else {
                            nuclide_string << n_vector_[j][2 * i] * Avogadro_Constant * lambda[i] * TOXCI[i] << " ";
                        }
                    }
                    if(n_vector_[size_tot][2 * i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][2 * i] * Avogadro_Constant * lambda[i] * TOXCI[i];
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    toxicity_node->InsertEndChild(nuclide);

                    /// 外回路
                    stringstream nuclide_loop_string;
                    nuclide_loop_string << scientific << setprecision(output_precision_);

                    for ( int j = 0; j < size_tot; ++j) {
                        if(n_vector_[j][2 * i + 1] < cutoff) {
                            nuclide_loop_string << 0.0 <<" ";
                        } else {
                            nuclide_loop_string << n_vector_[j][2 * i + 1] * Avogadro_Constant * lambda[i] * TOXCI[i] << " ";
                        }
                    }
                    if (n_vector_[size_tot][2 * i + 1] < cutoff) {
                        nuclide_loop_string << 0.0;
                    } else {
                        nuclide_loop_string << n_vector_[size_tot][2 * i + 1] * Avogadro_Constant * lambda[i] * TOXCI[i];
                    }

                    string nuclide_loop_s(nuclide_loop_string.str());
                    XMLText *nuclide_loop_text = doc.NewText(nuclide_loop_s.c_str());
                    nuclide_loop->InsertEndChild(nuclide_loop_text);
                    toxicity_loop_node->InsertEndChild(nuclide_loop);
                }
            }

        }

        if (if_print_fission_rate_ == true) {
            XMLElement* prod_node = doc.NewElement("NeuProdRate");
            prod_node->SetAttribute("unit","s-1");
            nuclides->InsertEndChild(prod_node);

            XMLElement *prod_tot = doc.NewElement("total");
            if (print_mode_ == 0) {
                stringstream prod_string;
                prod_string << scientific << setprecision(output_precision_);
                prod_string << prod_neu[0];
                prod_string << " ";
                prod_string << prod_neu.back();

                string prod_s(prod_string.str());

                XMLText *prod_text = doc.NewText(prod_s.c_str());
                prod_tot->InsertEndChild(prod_text);
                prod_node->InsertEndChild(prod_tot);
            } else {
                stringstream prod_string;
                prod_string << scientific << setprecision(output_precision_);
                for (int j = 0; j < size_tot; ++j) {
                    prod_string << prod_neu[j] << " ";
                }
                prod_string << prod_neu.back();

                string prod_s(prod_string.str());
                XMLText *prod_text = doc.NewText(prod_s.c_str());
                prod_tot->InsertEndChild(prod_text);
                prod_node->InsertEndChild(prod_tot);
            }

            for (int i = 0; i < nucl_size; ++ i) {
                XMLElement *nuclide = doc.NewElement("nuclide");
                int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                int NZ = nucl_id / 10000;
                int NA = (nucl_id - NZ * 10000) / 10;
                int NG = nucl_id - nucl_id / 10 * 10;
                string nucl_name;
                nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
                if(NG == 1) {
                    nucl_name += "m";
                }

                nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
                nuclide->SetAttribute("name",nucl_name.c_str());
                if(print_mode_ == 0) {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    nuclide_string << n_vector_[0][2 * i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24 * flux_vector_[0];
                    nuclide_string << " ";
                    if (n_vector_[size][2 * i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][2 * i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24 * flux_vector_[size_tot];
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    prod_node->InsertEndChild(nuclide);
                } else {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    for ( int j = 0; j < size_tot; ++j) {
                        if(n_vector_[j][2 * i] < cutoff) {
                            nuclide_string << 0.0 <<" ";
                        } else {
                            nuclide_string << n_vector_[j][2 * i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24 * flux_vector_[j] << " ";
                        }
                    }
                    if (n_vector_[size_tot][2 * i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][2 * i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24 * flux_vector_[size_tot];
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    prod_node->InsertEndChild(nuclide);
                }
            }

        }


        if( if_print_absorption_rate_ == true) {
            XMLElement* abs_node = doc.NewElement("NeuAbsRate");
            abs_node->SetAttribute("unit","s-1");
            nuclides->InsertEndChild(abs_node);

            XMLElement *abs_tot = doc.NewElement("total");
            if(print_mode_ == 0) {
                stringstream abs_string;
                abs_string << scientific << setprecision(output_precision_);
                abs_string << abs_neu[0];
                abs_string << " ";
                abs_string << abs_neu.back();

                string abs_s(abs_string.str());

                XMLText * abs_text = doc.NewText(abs_s.c_str());
                abs_tot->InsertEndChild(abs_text);
                abs_node->InsertEndChild(abs_tot);
            } else {
                stringstream abs_string;
                abs_string << scientific << setprecision(output_precision_);
                for(int j = 0; j < size_tot; ++j) {
                    abs_string << abs_neu[j] << " ";
                }
                abs_string << abs_neu.back();

                string abs_s(abs_string.str());
                XMLText *abs_text = doc.NewText(abs_s.c_str());
                abs_tot->InsertEndChild(abs_text);
                abs_node->InsertEndChild(abs_tot);
            }

            for (int i = 0; i < nucl_size; ++ i) {
                XMLElement *nuclide = doc.NewElement("nuclide");
                int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                int NZ = nucl_id / 10000;
                int NA = (nucl_id - NZ * 10000) / 10;
                int NG = nucl_id - nucl_id / 10 * 10;
                string nucl_name;
                nucl_name = ModecNuclideLibrary.element_name_list_[NZ - 1] + to_string(NA);
                if(NG == 1) {
                    nucl_name += "m";
                }

                nuclide->SetAttribute("zai",to_string(nucl_id).c_str());
                nuclide->SetAttribute("name",nucl_name.c_str());
                if(print_mode_ == 0) {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    nuclide_string << n_vector_[0][2 * i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24 * flux_vector_[0];
                    nuclide_string << " ";
                    if (n_vector_[size][2 * i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][2 * i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24 * flux_vector_[size_tot];
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText * nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    abs_node->InsertEndChild(nuclide);
                } else {
                    stringstream nuclide_string;
                    nuclide_string << scientific << setprecision(output_precision_);

                    for ( int j = 0; j < size_tot; ++j) {
                        if(n_vector_[j][2 * i] < cutoff) {
                            nuclide_string << 0.0 <<" ";
                        } else {
                            nuclide_string << n_vector_[j][2 * i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24 * flux_vector_[j] << " ";
                        }
                    }
                    if (n_vector_[size_tot][2 * i] < cutoff) {
                        nuclide_string << 0.0;
                    } else {
                        nuclide_string << n_vector_[size_tot][2 * i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24 * flux_vector_[size_tot];
                    }

                    string nuclide_s(nuclide_string.str());
                    XMLText *nuclide_text = doc.NewText(nuclide_s.c_str());
                    nuclide->InsertEndChild(nuclide_text);
                    abs_node->InsertEndChild(nuclide);
                }
            }
        }
    }

    if(input_filename_ == "") {
        output_filename_ =  input_filename_.substr(0, input_filename_.find_last_of(".")) + "_out" + ".xml";
    }
    doc.SaveFile(output_filename_.c_str());
}

void ModecClass::ModecOutput() {
    int Nucl_size = ModecNuclideLibrary.nuclide_number_;

    if( if_flow_mode_ == 0 ) { // 静止情况下的结果输出
        if (if_calculate_equilibrium_ == false) {
            if (input_filename_ != "modec.inp") {
                output_filename_ = input_filename_ + ".concentration";
            }

            output_file_ = work_direc_ + output_filename_;
            modec_out_.open(output_file_);
            modec_out_ << "*********************************** Nuclide Concentrations, " << dens_unit_ << " ***********************************";
            modec_out_ << '\n';
            modec_out_ << '\n';
            modec_out_.setf(ios::left);
            int size = evolution_mode_.size();
            int size_tot = 0;
            vector<string> time_unit;
            for (int i = 0; i < size; ++i) {
                size_tot += substep_[i];
            }

            vector<double> Cum_time;
            Cum_time.resize(size_tot + 1);

            int count = 0;
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < substep_[i]; ++j) {
                    Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                    time_unit.push_back(time_unit_[i]);
                    count++;
                }
            }


            modec_out_.width(8);
            modec_out_ << "";
            modec_out_.width(8);
            modec_out_ << "Time: ";
            modec_out_.width(18);
            modec_out_ << "Initial";

            if (print_mode_ == 0) {
                modec_out_.width(10);
                modec_out_.setf(ios::right);
                if (time_unit[size_tot - 1] == "d" || time_unit[size_tot - 1] == "day" || time_unit[size_tot - 1] == "days") {
                    Cum_time[size_tot] = Cum_time[size_tot] / 24.0 / 3600.0;
                }
                if (time_unit[size_tot - 1] == "m" || time_unit[size_tot - 1] == "month" || time_unit[size_tot - 1] == "months") {
                    Cum_time[size_tot] = Cum_time[size_tot] / 24.0 / 3600.0 / 30.0;
                } else if (time_unit[size_tot - 1] == "y" || time_unit[size_tot - 1] == "year" || time_unit[size_tot - 1] == "years") {
                    Cum_time[size_tot] = Cum_time[size_tot] / 24.0 / 3600.0 / 365.25;
                    time_unit[size_tot - 1] = "yr";
                }

                if (Cum_time[size_tot] >= 1.0e+3) {
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                } else {
                    modec_out_.precision(5);
                }

                modec_out_ << Cum_time[size_tot];
                modec_out_.unsetf(ios::right);
                modec_out_ << " ";
                modec_out_.width(7);
                modec_out_ << time_unit[size_tot - 1];
            } else {
                for (int j = 1; j <= size_tot; ++j) {
                    modec_out_.setf(ios::right);
                    if (time_unit[j - 1] == "d" || time_unit[j - 1] == "day" || time_unit[j - 1] == "days") {
                        Cum_time[j] = Cum_time[j] / 24.0 / 3600.0;
                    }
                    if (time_unit[j - 1] == "m" || time_unit[j - 1] == "month" || time_unit[j - 1] == "months") {
                        Cum_time[j] = Cum_time[j] / 24.0 / 3600.0 / 30.0;
                    } else if (time_unit[j - 1] == "y" || time_unit[j - 1] == "year" || time_unit[j - 1] == "years") {
                        Cum_time[j] = Cum_time[j] / 24.0 / 3600.0 / 365.25;
                        time_unit[j - 1] = "yr";
                    }
                    modec_out_.width(10);
                    if (Cum_time[j] >= 1.0e+3) {
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                    } else {
                        modec_out_.precision(5);
                    }
                    modec_out_ << Cum_time[j];
                    modec_out_.unsetf(ios::right);
                    modec_out_ << " ";
                    modec_out_.width(7);
                    modec_out_ << time_unit[j - 1];
                }
            }

            modec_out_ << '\n';

            modec_out_.width(8);
            modec_out_ << "";
            modec_out_.width(8);
            modec_out_ << "Power: ";

            if (print_mode_ == 0) {
                modec_out_.width(10);
                modec_out_.setf(ios::scientific | ios::uppercase);
                modec_out_.precision(4);
                modec_out_ << power_vector_[0];
                modec_out_.width(8);
                modec_out_ << "MW";

                modec_out_.width(10);
                modec_out_.setf(ios::scientific | ios::uppercase);
                modec_out_.precision(4);
                modec_out_ << power_vector_[size_tot];
                modec_out_.width(8);
                modec_out_ << "MW";
            } else {
                /*modec_out_.width(10);
                modec_out_.setf(ios::scientific | ios::uppercase);
                modec_out_.precision(4);
                modec_out_ << 0.0;
                modec_out_.width(8);
                modec_out_ << "MW";*/

                for (int j = 0; j <= size_tot; ++j) {
                    modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << power_vector_[j];
                    modec_out_.width(8);
                    modec_out_ << " MW";
                }
            }

            modec_out_ << '\n';

            modec_out_.width(8);
            modec_out_ << "";
            modec_out_.width(8);
            modec_out_ << "Flux: ";

            if (print_mode_ == 0) {
                modec_out_.width(10);
                modec_out_.setf(ios::scientific | ios::uppercase);
                modec_out_.precision(4);
                modec_out_ << flux_vector_[0];
                modec_out_.width(8);
                modec_out_ << "";

                modec_out_.width(10);
                modec_out_.setf(ios::scientific | ios::uppercase);
                modec_out_.precision(4);
                modec_out_ << flux_vector_[size_tot];
                modec_out_.width(8);
                modec_out_ << "";
            } else {
                //modec_out_.width(10);
                //modec_out_.setf(ios::scientific | ios::uppercase);
                //modec_out_.precision(4);
                //modec_out_ << 0.0;
                //modec_out_.width(8);
                //modec_out_ << "";
                for (int j = 0; j <= size_tot; ++j) {
                    modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << flux_vector_[j];
                    modec_out_.width(8);
                    modec_out_ << "";
                }
            }
            modec_out_ << '\n';

            //////////////    Calculate kinf_vector_, Neutron_Production_Rate and Neutron_Absorption_Rate      /////////
            vector<double> kinf;
            vector<double> prod_neu;
            vector<double> absorption_neu;
            if (if_print_kinf_ != 0 || if_print_fission_rate_ != 0 || if_print_absorption_rate_ != 0) {
                kinf.resize(size_tot + 1);
                prod_neu.resize(size_tot + 1);
                absorption_neu.resize(size_tot + 1);

                for (int j = 0; j <= size_tot; ++j) {
                    for (int i = 0; i < Nucl_size; ++i) {
                        prod_neu[j] += n_vector_[j][i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24;
                        absorption_neu[j] += n_vector_[j][i]*Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24;
                    }
                    prod_neu[j] = prod_neu[j] * flux_vector_[j];
                    absorption_neu[j] = absorption_neu[j] * flux_vector_[j];

                    if (absorption_neu[j] != 0.0) {
                        kinf[j] = prod_neu[j] / absorption_neu[j];
                    }

                }

                if (if_print_kinf_ == true) {

                    modec_out_.width(8);
                    modec_out_ << "";
                    modec_out_.width(8);
                    modec_out_ << "kinf: ";

                    if (print_mode_ == 0) {
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << kinf[0];
                        modec_out_.width(18);
                        modec_out_ << kinf[size_tot];
                        modec_out_ << '\n';
                    }
                    if (print_mode_ == 1) {
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            modec_out_ << kinf[j];
                        }
                        modec_out_ << '\n';
                    }
                }

                if (if_print_fission_rate_ != 0) {

                    modec_out_.width(2);
                    modec_out_ << "";
                    modec_out_.width(14);
                    modec_out_ << "NeuProdRate: ";

                    if (print_mode_ == 0) {
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << prod_neu[0];
                        modec_out_.width(18);
                        modec_out_ << prod_neu[size_tot];
                        modec_out_ << '\n';
                    }
                    if (print_mode_ == 1) {
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            modec_out_ << prod_neu[j];
                        }
                        modec_out_ << '\n';
                    }
                }

                if (if_print_absorption_rate_ != 0) {

                    modec_out_.width(2);
                    modec_out_ << "";
                    modec_out_.width(14);
                    modec_out_ << "NeuAbsRate: ";

                    if (print_mode_ == 0) {
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << absorption_neu[0];
                        modec_out_.width(18);
                        modec_out_ << absorption_neu[size_tot];
                        modec_out_ << '\n';
                    }
                    if (print_mode_ == 1) {
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            modec_out_ << absorption_neu[j];
                        }
                        modec_out_ << '\n';
                    }
                }
                modec_out_ << '\n';
            }
            //////////////////////////////////////////////////////////////////////////////////////////////////


            modec_out_ << '\n';

            for (unsigned int i = 0; i < Nucl_size; ++i) {
                modec_out_.width(8);
                int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                modec_out_ << nucl_id;
                int NZ = nucl_id / 10000;
                int NA = (nucl_id - NZ * 10000) / 10;
                stringstream ss;
                ss << NA;
                int NG = nucl_id - nucl_id / 10 * 10;
                string Nucl_Name;
                Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                if (NG == 1) {
                    Nucl_Name += "m";
                }
                modec_out_.width(8);
                modec_out_ << Nucl_Name;

                double convert_coeff;

                if (dens_unit_ == "mol") {
                    convert_coeff = 1.0;
                } else if (dens_unit_ == "g") {
                    convert_coeff = double(NA);
                } else if (dens_unit_ == "kg") {
                    convert_coeff = double(NA) / 1000.0;
                } else if (dens_unit_ == "atom") {
                    convert_coeff = Avogadro_Constant;
                } else if (dens_unit_ == "atom/(barn-cm)") {
                    convert_coeff = Avogadro_Constant * (1.0e-24);
                }

                if (print_mode_ == 0) {
                    modec_out_.width(18);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(9);
                    modec_out_ << n_vector_[0][i] * convert_coeff;
                    modec_out_.width(18);
                    if (n_vector_[size_tot][i] < cutoff) {
                        modec_out_ << 0.0;
                    } else {
                        modec_out_ << n_vector_[size_tot][i] * convert_coeff;
                    }

                    modec_out_ << '\n';
                }
                if (print_mode_ == 1) {
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(9);
                    for (int j = 0; j <= size_tot; ++j) {
                        modec_out_.width(18);
                        if (n_vector_[j][i] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[j][i] * convert_coeff;
                        }
                    }
                    modec_out_ << '\n';
                }
            }

            modec_out_.close();

            if (if_print_activity_ == true) {

                vector<double > lamda;
                lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

                output_filename_ = input_filename_ + ".radioactivity";

                output_file_ = work_direc_ + output_filename_;
                modec_out_.open(output_file_);
                modec_out_ << "*********************************** Nuclide Radioactivity, Ci ***********************************";
                modec_out_ << '\n';
                modec_out_ << '\n';

                modec_out_.setf(ios::left);
                int size = evolution_mode_.size();
                int size_tot = 0;
                for (int i = 0; i < size; ++i) {
                    size_tot += substep_[i];
                }

                vector<double> Cum_time;
                Cum_time.resize(size_tot + 1);

                int count = 0;
                for (int i = 0; i < size; ++i) {
                    for (int j = 0; j < substep_[i]; ++j) {
                        Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                        count++;
                    }
                }


                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(18);
                modec_out_ << "Initial";

                if (print_mode_ == 0) {
                    modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << Cum_time[size_tot];
                    modec_out_.width(8);
                    modec_out_ << "s";
                } else {
                    for (int j = 1; j <= size_tot; ++j) {
                        modec_out_.width(10);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                        modec_out_ << Cum_time[j];
                        modec_out_.width(8);
                        modec_out_ << "s";
                    }
                }

                modec_out_ << '\n';
                modec_out_ << '\n';

                for (unsigned int i = 0; i < Nucl_size; ++i) {
                    modec_out_.width(8);
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    modec_out_ << nucl_id;
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    stringstream ss;
                    ss << NA;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string Nucl_Name;
                    Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                    if (NG == 1) {
                        Nucl_Name += "m";
                    }
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;
                    if (print_mode_ == 0) {
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << n_vector_[0][i]*Avogadro_Constant*lamda[i] / (3.7e10);
                        modec_out_.width(18);
                        if (n_vector_[size_tot][i] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[size_tot][i]*Avogadro_Constant*lamda[i] / (3.7e10); // 单位为Ci
                        }

                        modec_out_ << '\n';
                    }
                    if (print_mode_ == 1) {
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            if (n_vector_[j][i] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[j][i]*Avogadro_Constant*lamda[i] / (3.7e10); // 单位为Ci
                            }
                        }
                        modec_out_ << '\n';
                    }
                }

                modec_out_.close();

            }

            if (if_print_decayenergy_ == true) {
                vector<double > lamda;
                lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

                vector<double > Q;
                Q = ModecNuclideLibrary.nuclide_library_vector_[3];

                output_filename_ = input_filename_ + ".decayheat";

                output_file_ = work_direc_ + output_filename_;
                modec_out_.open(output_file_);
                modec_out_ << "*********************************** Nuclide Decay Heat, Watts ***********************************";
                modec_out_ << '\n';
                modec_out_ << '\n';

                modec_out_.setf(ios::left);
                int size = evolution_mode_.size();
                int size_tot = 0;
                for (int i = 0; i < size; ++i) {
                    size_tot += substep_[i];
                }

                vector<double> Cum_time;
                Cum_time.resize(size_tot + 1);

                int count = 0;
                for (int i = 0; i < size; ++i) {
                    for (int j = 0; j < substep_[i]; ++j) {
                        Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                        count++;
                    }
                }


                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(18);
                modec_out_ << "Initial";

                if (print_mode_ == 0) {
                    modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << Cum_time[size_tot];
                    modec_out_.width(8);
                    modec_out_ << "s";
                } else {
                    for (int j = 1; j <= size_tot; ++j) {
                        modec_out_.width(10);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                        modec_out_ << Cum_time[j];
                        modec_out_.width(8);
                        modec_out_ << "s";
                    }
                }

                modec_out_ << '\n';
                modec_out_ << '\n';

                for (unsigned int i = 0; i < Nucl_size; ++i) {
                    modec_out_.width(8);
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    modec_out_ << nucl_id;
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    stringstream ss;
                    ss << NA;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string Nucl_Name;
                    Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                    if (NG == 1) {
                        Nucl_Name += "m";
                    }
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;
                    if (print_mode_ == 0) {
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << n_vector_[0][i]*Avogadro_Constant*lamda[i] * Q[i] * Electron_Coulomb*1.0e6;
                        modec_out_.width(18);
                        if (n_vector_[size_tot][i] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[size_tot][i]*Avogadro_Constant*lamda[i] * Q[i] * Electron_Coulomb*1.0e6; // 转换成单位 W
                        }

                        modec_out_ << '\n';
                    }
                    if (print_mode_ == 1) {
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            if (n_vector_[j][i] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[j][i]*Avogadro_Constant*lamda[i]*Q[i]* Electron_Coulomb*1.0e6;
                            }
                        }
                        modec_out_ << '\n';
                    }
                }

                modec_out_.close();

            }

            if (if_print_ampc_ == true) {
                vector<double > lamda;
                lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

                vector<double > AMPC;
                AMPC = ModecNuclideLibrary.nuclide_library_vector_[4];

                output_filename_ = input_filename_ + ".ampc";

                output_file_ = work_direc_ + output_filename_;
                modec_out_.open(output_file_);
                modec_out_ << "*********************************** Nuclide Radioactivity, m3-air/Bq ***********************************";
                modec_out_ << '\n';
                modec_out_ << '\n';

                modec_out_.setf(ios::left);
                int size = evolution_mode_.size();
                int size_tot = 0;
                for (int i = 0; i < size; ++i) {
                    size_tot += substep_[i];
                }

                vector<double> Cum_time;
                Cum_time.resize(size_tot + 1);

                int count = 0;
                for (int i = 0; i < size; ++i) {
                    for (int j = 0; j < substep_[i]; ++j) {
                        Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                        count++;
                    }
                }


                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(18);
                modec_out_ << "Initial";

                if (print_mode_ == 0) {
                    modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << Cum_time[size_tot];
                    modec_out_.width(8);
                    modec_out_ << "s";
                } else {
                    for (int j = 1; j <= size_tot; ++j) {
                        modec_out_.width(10);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                        modec_out_ << Cum_time[j];
                        modec_out_.width(8);
                        modec_out_ << "s";
                    }
                }

                modec_out_ << '\n';
                modec_out_ << '\n';

                for (unsigned int i = 0; i < Nucl_size; ++i) {
                    modec_out_.width(8);
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    modec_out_ << nucl_id;
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    stringstream ss;
                    ss << NA;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string Nucl_Name;
                    Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                    if (NG == 1) {
                        Nucl_Name += "m";
                    }
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;
                    if (print_mode_ == 0) {
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << n_vector_[0][i]*Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
                        modec_out_.width(18);
                        if (n_vector_[size_tot][i] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[size_tot][i]*Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
                        }

                        modec_out_ << '\n';
                    }
                    if (print_mode_ == 1) {
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            if (n_vector_[j][i] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[j][i]*Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
                            }
                        }
                        modec_out_ << '\n';
                    }
                }

                modec_out_.close();

            }

            if (if_print_wmpc_ == true) {
                vector<double > lamda;
                lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

                vector<double > WMPC;
                WMPC = ModecNuclideLibrary.nuclide_library_vector_[4];

                output_filename_ = input_filename_ + ".wmpc";

                output_file_ = work_direc_ + output_filename_;
                modec_out_.open(output_file_);
                modec_out_ << "*********************************** Nuclide Radioactivity, m3-air/Bq ***********************************";
                modec_out_ << '\n';
                modec_out_ << '\n';

                modec_out_.setf(ios::left);
                int size = evolution_mode_.size();
                int size_tot = 0;
                for (int i = 0; i < size; ++i) {
                    size_tot += substep_[i];
                }

                vector<double> Cum_time;
                Cum_time.resize(size_tot + 1);

                int count = 0;
                for (int i = 0; i < size; ++i) {
                    for (int j = 0; j < substep_[i]; ++j) {
                        Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                        count++;
                    }
                }


                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(18);
                modec_out_ << "Initial";

                if (print_mode_ == 0) {
                    modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << Cum_time[size_tot];
                    modec_out_.width(8);
                    modec_out_ << "s";
                } else {
                    for (int j = 1; j <= size_tot; ++j) {
                        modec_out_.width(10);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                        modec_out_ << Cum_time[j];
                        modec_out_.width(8);
                        modec_out_ << "s";
                    }
                }

                modec_out_ << '\n';
                modec_out_ << '\n';

                for (unsigned int i = 0; i < Nucl_size; ++i) {
                    modec_out_.width(8);
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    modec_out_ << nucl_id;
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    stringstream ss;
                    ss << NA;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string Nucl_Name;
                    Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                    if (NG == 1) {
                        Nucl_Name += "m";
                    }
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;
                    if (print_mode_ == 0) {
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << n_vector_[0][i]*Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
                        modec_out_.width(18);
                        if (n_vector_[size_tot][i] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[size_tot][i]*Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
                        }

                        modec_out_ << '\n';
                    }
                    if (print_mode_ == 1) {
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            if (n_vector_[j][i] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[j][i]*Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
                            }
                        }
                        modec_out_ << '\n';
                    }
                }

                modec_out_.close();

            }

            if (if_print_toxicity_ == true) {
                vector<double > lamda;
                lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

                vector<double > TOXCI;
                TOXCI = ModecNuclideLibrary.nuclide_library_vector_[10];

                output_filename_ = input_filename_ + ".toxicity";

                output_file_ = work_direc_ + output_filename_;
                modec_out_.open(output_file_);
                modec_out_ << "*********************************** Nuclide Radiotoxicity, Sv ***********************************";
                modec_out_ << '\n';
                modec_out_ << '\n';

                modec_out_.setf(ios::left);
                int size = evolution_mode_.size();
                int size_tot = 0;
                for (int i = 0; i < size; ++i) {
                    size_tot += substep_[i];
                }

                vector<double> Cum_time;
                Cum_time.resize(size_tot + 1);

                int count = 0;
                for (int i = 0; i < size; ++i) {
                    for (int j = 0; j < substep_[i]; ++j) {
                        Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                        count++;
                    }
                }


                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(18);
                modec_out_ << "Initial";

                if (print_mode_ == 0) {
                    modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << Cum_time[size_tot];
                    modec_out_.width(8);
                    modec_out_ << "s";
                } else {
                    for (int j = 1; j <= size_tot; ++j) {
                        modec_out_.width(10);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                        modec_out_ << Cum_time[j];
                        modec_out_.width(8);
                        modec_out_ << "s";
                    }
                }

                modec_out_ << '\n';
                modec_out_ << '\n';

                for (unsigned int i = 0; i < Nucl_size; ++i) {
                    modec_out_.width(8);
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    modec_out_ << nucl_id;
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    stringstream ss;
                    ss << NA;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string Nucl_Name;
                    Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                    if (NG == 1) {
                        Nucl_Name += "m";
                    }
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;
                    if (print_mode_ == 0) {
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << n_vector_[0][i]*Avogadro_Constant*lamda[i] * TOXCI[i];
                        modec_out_.width(18);
                        if (n_vector_[size_tot][i] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[size_tot][i]*Avogadro_Constant*lamda[i] * TOXCI[i];
                        }

                        modec_out_ << '\n';
                    }
                    if (print_mode_ == 1) {
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            if (n_vector_[j][i] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[j][i]*Avogadro_Constant*lamda[i] * TOXCI[i];
                            }
                        }
                        modec_out_ << '\n';
                    }
                }

                modec_out_.close();

            }

            if (if_print_fission_rate_ == 2) {
                output_filename_ = input_filename_ + ".NeuProdRate";

                output_file_ = work_direc_ + output_filename_;
                modec_out_.open(output_file_);
                modec_out_ << "*********************************** Neutron Production Rate, s-1***********************************";
                modec_out_ << '\n';
                modec_out_ << '\n';

                modec_out_.setf(ios::left);
                int size = evolution_mode_.size();
                int size_tot = 0;
                for (int i = 0; i < size; ++i) {
                    size_tot += substep_[i];
                }

                vector<double> Cum_time;
                Cum_time.resize(size_tot + 1);

                int count = 0;
                for (int i = 0; i < size; ++i) {
                    for (int j = 0; j < substep_[i]; ++j) {
                        Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                        count++;
                    }
                }


                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(18);
                modec_out_ << "Initial";

                if (print_mode_ == 0) {
                    modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << Cum_time[size_tot];
                    modec_out_.width(8);
                    modec_out_ << "s";
                } else {
                    for (int j = 1; j <= size_tot; ++j) {
                        modec_out_.width(10);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                        modec_out_ << Cum_time[j];
                        modec_out_.width(8);
                        modec_out_ << "s";
                    }
                }

                modec_out_ << '\n';
                modec_out_ << '\n';

                for (unsigned int i = 0; i < Nucl_size; ++i) {
                    modec_out_.width(8);
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    modec_out_ << nucl_id;
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    stringstream ss;
                    ss << NA;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string Nucl_Name;
                    Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                    if (NG == 1) {
                        Nucl_Name += "m";
                    }
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;
                    if (print_mode_ == 0) {
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << n_vector_[0][i]*Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24 * flux_vector_[0];
                        modec_out_.width(18);

                        if (n_vector_[size_tot][i] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[size_tot][i]*Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24* flux_vector_[size_tot];
                        }
                        modec_out_ << '\n';
                    }
                    if (print_mode_ == 1) {
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            if (n_vector_[j][i] < cutoff || n_vector_[j][i] < 0.0) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[j][i]*Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24* flux_vector_[j]; // 单位为Ci
                            }
                        }
                        modec_out_ << '\n';
                    }
                }
                modec_out_ << '\n';
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "Total: ";
                for (int j = 0; j <= size_tot; ++j) {
                    modec_out_.width(18);

                    modec_out_ << prod_neu[j]; // 单位为Ci

                }
                modec_out_ << '\n';

                modec_out_.close();

            }

            if (if_print_absorption_rate_ == 2) {
                output_filename_ = input_filename_ + ".NeuAbsRate";

                output_file_ = work_direc_ + output_filename_;
                modec_out_.open(output_file_);
                modec_out_ << "*********************************** Neutron Absorption Rate, s-1 ***********************************";
                modec_out_ << '\n';
                modec_out_ << '\n';

                modec_out_.setf(ios::left);
                int size = evolution_mode_.size();
                int size_tot = 0;
                for (int i = 0; i < size; ++i) {
                    size_tot += substep_[i];
                }

                vector<double> Cum_time;
                Cum_time.resize(size_tot + 1);

                int count = 0;
                for (int i = 0; i < size; ++i) {
                    for (int j = 0; j < substep_[i]; ++j) {
                        Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                        count++;
                    }
                }


                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(18);
                modec_out_ << "Initial";

                if (print_mode_ == 0) {
                    modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << Cum_time[size_tot];
                    modec_out_.width(8);
                    modec_out_ << "s";
                } else {
                    for (int j = 1; j <= size_tot; ++j) {
                        modec_out_.width(10);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                        modec_out_ << Cum_time[j];
                        modec_out_.width(8);
                        modec_out_ << "s";
                    }
                }

                modec_out_ << '\n';
                modec_out_ << '\n';

                for (unsigned int i = 0; i < Nucl_size; ++i) {
                    modec_out_.width(8);
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    modec_out_ << nucl_id;
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    stringstream ss;
                    ss << NA;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string Nucl_Name;
                    Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                    if (NG == 1) {
                        Nucl_Name += "m";
                    }
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;
                    if (print_mode_ == 0) {
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << n_vector_[0][i]*Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24* flux_vector_[0];
                        modec_out_.width(18);
                        if (n_vector_[size_tot][i] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[size_tot][i]*Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24* flux_vector_[size_tot];
                        }
                        modec_out_ << '\n';
                    }
                    if (print_mode_ == 1) {
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            if (n_vector_[j][i] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[j][i]*Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24* flux_vector_[j];
                            }
                        }
                        modec_out_ << '\n';
                    }
                }

                modec_out_ << '\n';
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "Total: ";
                for (int j = 0; j <= size_tot; ++j) {
                    modec_out_.width(18);

                    modec_out_ << absorption_neu[j]; // 单位为Ci

                }
                modec_out_ << '\n';

                modec_out_.close();

            }
            if (if_tracking_stockage == true) {
                output_filename_ = input_filename_ + ".stockage.concentration";

                output_file_ = work_direc_ + output_filename_;
                modec_out_.open(output_file_);
                modec_out_ << "*********************************** Nuclide Concentrations, mol ***********************************";
                modec_out_ << '\n';
                modec_out_ << '\n';
                modec_out_.setf(ios::left);
                int size = evolution_mode_.size();
                int size_tot = 0;
                vector<string> time_unit;
                for (int i = 0; i < size; ++i) {
                    size_tot += substep_[i];
                }

                vector<double> Cum_time;
                Cum_time.resize(size_tot + 1);

                int count = 0;
                for (int i = 0; i < size; ++i) {
                    for (int j = 0; j < substep_[i]; ++j) {
                        Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                        time_unit.push_back(time_unit_[i]);
                        count++;
                    }
                }


                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "Time: ";
                modec_out_.width(18);
                modec_out_ << "Initial";

                if (print_mode_ == 0) {
                    modec_out_.width(10);
                    modec_out_.setf(ios::right);
                    if (time_unit[size_tot - 1] == "d" || time_unit[size_tot - 1] == "day" || time_unit[size_tot - 1] == "days") {
                        Cum_time[size_tot] = Cum_time[size_tot] / 24.0 / 3600.0;
                    }
                    if (time_unit[size_tot - 1] == "m" || time_unit[size_tot - 1] == "month" || time_unit[size_tot - 1] == "months") {
                        Cum_time[size_tot] = Cum_time[size_tot] / 24.0 / 3600.0 / 30.0;
                    } else if (time_unit[size_tot - 1] == "y" || time_unit[size_tot - 1] == "year" || time_unit[size_tot - 1] == "years") {
                        Cum_time[size_tot] = Cum_time[size_tot] / 24.0 / 3600.0 / 365.25;
                        time_unit[size_tot - 1] = "yr";
                    }

                    if (Cum_time[size_tot] >= 1.0e+3) {
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                    } else {
                        modec_out_.unsetf(ios::scientific | ios::uppercase);
                        modec_out_.precision(5);
                    }

                    modec_out_ << Cum_time[size_tot];
                    modec_out_.unsetf(ios::right);
                    modec_out_ << " ";
                    modec_out_.width(7);
                    modec_out_ << time_unit[size_tot - 1];
                } else {
                    for (int j = 1; j <= size_tot; ++j) {
                        modec_out_.setf(ios::right);
                        if (time_unit[j - 1] == "d" || time_unit[j - 1] == "day" || time_unit[j - 1] == "days") {
                            Cum_time[j] = Cum_time[j] / 24.0 / 3600.0;
                        }
                        if (time_unit[j - 1] == "m" || time_unit[j - 1] == "month" || time_unit[j - 1] == "months") {
                            Cum_time[j] = Cum_time[j] / 24.0 / 3600.0 / 30.0;
                        } else if (time_unit[j - 1] == "y" || time_unit[j - 1] == "year" || time_unit[j - 1] == "years") {
                            Cum_time[j] = Cum_time[j] / 24.0 / 3600.0 / 365.25;
                            time_unit[j - 1] = "yr";
                        }
                        modec_out_.width(10);
                        if (Cum_time[j] >= 1.0e+3) {
                            modec_out_.setf(ios::scientific | ios::uppercase);
                            modec_out_.precision(4);
                        } else {
                            modec_out_.unsetf(ios::scientific | ios::uppercase);
                            modec_out_.precision(5);
                        }

                        modec_out_ << Cum_time[j];
                        modec_out_.unsetf(ios::right);
                        modec_out_ << " ";
                        modec_out_.width(7);
                        modec_out_ << time_unit[j - 1];
                    }
                }

                modec_out_ << '\n';

                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "Power: ";

                if (print_mode_ == 0) {
                    modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << 0.0;
                    modec_out_.width(8);
                    modec_out_ << "MW";

                    modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << 0.0;
                    modec_out_.width(8);
                    modec_out_ << "MW";
                } else {
                    /*modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << 0.0;
                    modec_out_.width(8);
                    modec_out_ << "MW";*/

                    for (int j = 0; j <= size_tot; ++j) {
                        modec_out_.width(10);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                        modec_out_ << 0.0;
                        modec_out_.width(8);
                        modec_out_ << " MW";
                    }
                }

                modec_out_ << '\n';

                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "flux_: ";

                if (print_mode_ == 0) {
                    modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << 0.0;
                    modec_out_.width(8);
                    modec_out_ << "";

                    modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << 0.0;
                    modec_out_.width(8);
                    modec_out_ << "";
                } else {
                    //modec_out_.width(10);
                    //modec_out_.setf(ios::scientific | ios::uppercase);
                    //modec_out_.precision(4);
                    //modec_out_ << 0.0;
                    //modec_out_.width(8);
                    //modec_out_ << "";
                    for (int j = 0; j <= size_tot; ++j) {
                        modec_out_.width(10);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                        modec_out_ << 0.0;
                        modec_out_.width(8);
                        modec_out_ << "";
                    }
                }
                modec_out_ << '\n';
                modec_out_ << '\n';

                for (unsigned int i = 0; i < Nucl_size; ++i) {
                    modec_out_.width(8);
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    modec_out_ << nucl_id;
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    stringstream ss;
                    ss << NA;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string Nucl_Name;
                    Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                    if (NG == 1) {
                        Nucl_Name += "m";
                    }
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;
                    if (print_mode_ == 0) {
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << n_vector_[0][i + Nucl_size];
                        modec_out_.width(18);
                        if (n_vector_[size_tot][i + Nucl_size] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[size_tot][i + Nucl_size];
                        }

                        modec_out_ << '\n';
                    }
                    if (print_mode_ == 1) {
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            if (n_vector_[j][i + Nucl_size] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[j][i + Nucl_size];
                            }
                        }
                        modec_out_ << '\n';
                    }
                }

                modec_out_.close();

                if (if_print_stockage_activity_ == true) {

                    vector<double > lamda;
                    lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

                    output_filename_ = input_filename_ + ".stockage.radioactivity";

                    output_file_ = work_direc_ + output_filename_;
                    modec_out_.open(output_file_);
                    modec_out_ << "*********************************** Nuclide Radioactivity, Ci ***********************************";
                    modec_out_ << '\n';
                    modec_out_ << '\n';

                    modec_out_.setf(ios::left);
                    int size = evolution_mode_.size();
                    int size_tot = 0;
                    for (int i = 0; i < size; ++i) {
                        size_tot += substep_[i];
                    }

                    vector<double> Cum_time;
                    Cum_time.resize(size_tot + 1);

                    int count = 0;
                    for (int i = 0; i < size; ++i) {
                        for (int j = 0; j < substep_[i]; ++j) {
                            Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                            count++;
                        }
                    }


                    modec_out_.width(8);
                    modec_out_ << "";
                    modec_out_.width(8);
                    modec_out_ << "";
                    modec_out_.width(18);
                    modec_out_ << "Initial";

                    if (print_mode_ == 0) {
                        modec_out_.width(10);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                        modec_out_ << Cum_time[size_tot];
                        modec_out_.width(8);
                        modec_out_ << "s";
                    } else {
                        for (int j = 1; j <= size_tot; ++j) {
                            modec_out_.width(10);
                            modec_out_.setf(ios::scientific | ios::uppercase);
                            modec_out_.precision(4);
                            modec_out_ << Cum_time[j];
                            modec_out_.width(8);
                            modec_out_ << "s";
                        }
                    }

                    modec_out_ << '\n';
                    modec_out_ << '\n';

                    for (unsigned int i = 0; i < Nucl_size; ++i) {
                        modec_out_.width(8);
                        int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                        modec_out_ << nucl_id;
                        int NZ = nucl_id / 10000;
                        int NA = (nucl_id - NZ * 10000) / 10;
                        stringstream ss;
                        ss << NA;
                        int NG = nucl_id - nucl_id / 10 * 10;
                        string Nucl_Name;
                        Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                        if (NG == 1) {
                            Nucl_Name += "m";
                        }
                        modec_out_.width(8);
                        modec_out_ << Nucl_Name;
                        if (print_mode_ == 0) {
                            modec_out_.width(18);
                            modec_out_.setf(ios::scientific | ios::uppercase);
                            modec_out_.precision(9);
                            modec_out_ << n_vector_[0][i + Nucl_size]*Avogadro_Constant*lamda[i] / (3.7e10);
                            modec_out_.width(18);
                            if (n_vector_[size_tot][i + Nucl_size] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[size_tot][i + Nucl_size]*Avogadro_Constant*lamda[i] / (3.7e10); // 单位为Ci
                            }

                            modec_out_ << '\n';
                        }
                        if (print_mode_ == 1) {
                            modec_out_.setf(ios::scientific | ios::uppercase);
                            modec_out_.precision(9);
                            for (int j = 0; j <= size_tot; ++j) {
                                modec_out_.width(18);
                                if (n_vector_[j][i + Nucl_size] < cutoff) {
                                    modec_out_ << 0.0;
                                } else {
                                    modec_out_ << n_vector_[j][i + Nucl_size]*Avogadro_Constant*lamda[i] / (3.7e10); // 单位为Ci
                                }
                            }
                            modec_out_ << '\n';
                        }
                    }

                    modec_out_.close();

                }

                if (if_print_stockage_decayenergy_ == true) {
                    vector<double > lamda;
                    lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

                    vector<double > Q;
                    Q = ModecNuclideLibrary.nuclide_library_vector_[3];

                    output_filename_ = input_filename_ + ".stockage.decayheat";
                    output_file_ = work_direc_ + output_filename_;
                    modec_out_.open(output_file_);
                    modec_out_ << "*********************************** Nuclide Decay Heat, Watts ***********************************";
                    modec_out_ << '\n';
                    modec_out_ << '\n';

                    modec_out_.setf(ios::left);
                    int size = evolution_mode_.size();
                    int size_tot = 0;
                    for (int i = 0; i < size; ++i) {
                        size_tot += substep_[i];
                    }

                    vector<double> Cum_time;
                    Cum_time.resize(size_tot + 1);

                    int count = 0;
                    for (int i = 0; i < size; ++i) {
                        for (int j = 0; j < substep_[i]; ++j) {
                            Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                            count++;
                        }
                    }


                    modec_out_.width(8);
                    modec_out_ << "";
                    modec_out_.width(8);
                    modec_out_ << "";
                    modec_out_.width(18);
                    modec_out_ << "Initial";

                    if (print_mode_ == 0) {
                        modec_out_.width(10);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                        modec_out_ << Cum_time[size_tot];
                        modec_out_.width(8);
                        modec_out_ << "s";
                    } else {
                        for (int j = 1; j <= size_tot; ++j) {
                            modec_out_.width(10);
                            modec_out_.setf(ios::scientific | ios::uppercase);
                            modec_out_.precision(4);
                            modec_out_ << Cum_time[j];
                            modec_out_.width(8);
                            modec_out_ << "s";
                        }
                    }

                    modec_out_ << '\n';
                    modec_out_ << '\n';

                    for (unsigned int i = 0; i < Nucl_size; ++i) {
                        modec_out_.width(8);
                        int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                        modec_out_ << nucl_id;
                        int NZ = nucl_id / 10000;
                        int NA = (nucl_id - NZ * 10000) / 10;
                        stringstream ss;
                        ss << NA;
                        int NG = nucl_id - nucl_id / 10 * 10;
                        string Nucl_Name;
                        Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                        if (NG == 1) {
                            Nucl_Name += "m";
                        }
                        modec_out_.width(8);
                        modec_out_ << Nucl_Name;
                        if (print_mode_ == 0) {
                            modec_out_.width(18);
                            modec_out_.setf(ios::scientific | ios::uppercase);
                            modec_out_.precision(9);
                            modec_out_ << n_vector_[0][i + Nucl_size]*Avogadro_Constant*lamda[i] * Q[i] * Electron_Coulomb*1.0e6;
                            modec_out_.width(18);
                            if (n_vector_[size_tot][i + Nucl_size] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[size_tot][i + Nucl_size]*Avogadro_Constant*lamda[i] * Q[i] * Electron_Coulomb*1.0e6; // 转换成单位 W
                            }

                            modec_out_ << '\n';
                        }
                        if (print_mode_ == 1) {
                            modec_out_.setf(ios::scientific | ios::uppercase);
                            modec_out_.precision(9);
                            for (int j = 0; j <= size_tot; ++j) {
                                modec_out_.width(18);
                                if (n_vector_[j][i + Nucl_size] < cutoff) {
                                    modec_out_ << 0.0;
                                } else {
                                    modec_out_ << n_vector_[j][i + Nucl_size]*Avogadro_Constant*lamda[i]*Q[i]* Electron_Coulomb*1.0e6;
                                }
                            }
                            modec_out_ << '\n';
                        }
                    }

                    modec_out_.close();

                }

                if (if_print_stockage_ampc_ == true) {
                    vector<double > lamda;
                    lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

                    vector<double > AMPC;
                    AMPC = ModecNuclideLibrary.nuclide_library_vector_[4];

                    output_filename_ = input_filename_ + ".stockage.ampc";
                    output_file_ = work_direc_ + output_filename_;
                    modec_out_.open(output_file_);
                    modec_out_ << "*********************************** Nuclide Radioactivity, m3-air/Bq ***********************************";
                    modec_out_ << '\n';
                    modec_out_ << '\n';

                    modec_out_.setf(ios::left);
                    int size = evolution_mode_.size();
                    int size_tot = 0;
                    for (int i = 0; i < size; ++i) {
                        size_tot += substep_[i];
                    }

                    vector<double> Cum_time;
                    Cum_time.resize(size_tot + 1);

                    int count = 0;
                    for (int i = 0; i < size; ++i) {
                        for (int j = 0; j < substep_[i]; ++j) {
                            Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                            count++;
                        }
                    }


                    modec_out_.width(8);
                    modec_out_ << "";
                    modec_out_.width(8);
                    modec_out_ << "";
                    modec_out_.width(18);
                    modec_out_ << "Initial";

                    if (print_mode_ == 0) {
                        modec_out_.width(10);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                        modec_out_ << Cum_time[size_tot];
                        modec_out_.width(8);
                        modec_out_ << "s";
                    } else {
                        for (int j = 1; j <= size_tot; ++j) {
                            modec_out_.width(10);
                            modec_out_.setf(ios::scientific | ios::uppercase);
                            modec_out_.precision(4);
                            modec_out_ << Cum_time[j];
                            modec_out_.width(8);
                            modec_out_ << "s";
                        }
                    }

                    modec_out_ << '\n';
                    modec_out_ << '\n';

                    for (unsigned int i = 0; i < Nucl_size; ++i) {
                        modec_out_.width(8);
                        int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                        modec_out_ << nucl_id;
                        int NZ = nucl_id / 10000;
                        int NA = (nucl_id - NZ * 10000) / 10;
                        stringstream ss;
                        ss << NA;
                        int NG = nucl_id - nucl_id / 10 * 10;
                        string Nucl_Name;
                        Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                        if (NG == 1) {
                            Nucl_Name += "m";
                        }
                        modec_out_.width(8);
                        modec_out_ << Nucl_Name;
                        if (print_mode_ == 0) {
                            modec_out_.width(18);
                            modec_out_.setf(ios::scientific | ios::uppercase);
                            modec_out_.precision(9);
                            modec_out_ << n_vector_[0][i + Nucl_size]*Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
                            modec_out_.width(18);
                            if (n_vector_[size_tot][i + Nucl_size] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[size_tot][i + Nucl_size]*Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
                            }

                            modec_out_ << '\n';
                        }
                        if (print_mode_ == 1) {
                            modec_out_.setf(ios::scientific | ios::uppercase);
                            modec_out_.precision(9);
                            for (int j = 0; j <= size_tot; ++j) {
                                modec_out_.width(18);
                                if (n_vector_[j][i + Nucl_size] < cutoff) {
                                    modec_out_ << 0.0;
                                } else {
                                    modec_out_ << n_vector_[j][i + Nucl_size]*Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
                                }
                            }
                            modec_out_ << '\n';
                        }
                    }

                    modec_out_.close();

                }

                if (if_print_stockage_wmpc_ == true) {
                    vector<double > lamda;
                    lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

                    vector<double > WMPC;
                    WMPC = ModecNuclideLibrary.nuclide_library_vector_[4];

                    output_filename_ = input_filename_ + ".stockage.wmpc";
                    output_file_ = work_direc_ + output_filename_;
                    modec_out_.open(output_file_);
                    modec_out_ << "*********************************** Nuclide Radioactivity, m3-air/Bq ***********************************";
                    modec_out_ << '\n';
                    modec_out_ << '\n';

                    modec_out_.setf(ios::left);
                    int size = evolution_mode_.size();
                    int size_tot = 0;
                    for (int i = 0; i < size; ++i) {
                        size_tot += substep_[i];
                    }

                    vector<double> Cum_time;
                    Cum_time.resize(size_tot + 1);

                    int count = 0;
                    for (int i = 0; i < size; ++i) {
                        for (int j = 0; j < substep_[i]; ++j) {
                            Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                            count++;
                        }
                    }


                    modec_out_.width(8);
                    modec_out_ << "";
                    modec_out_.width(8);
                    modec_out_ << "";
                    modec_out_.width(18);
                    modec_out_ << "Initial";

                    if (print_mode_ == 0) {
                        modec_out_.width(10);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                        modec_out_ << Cum_time[size_tot];
                        modec_out_.width(8);
                        modec_out_ << "s";
                    } else {
                        for (int j = 1; j <= size_tot; ++j) {
                            modec_out_.width(10);
                            modec_out_.setf(ios::scientific | ios::uppercase);
                            modec_out_.precision(4);
                            modec_out_ << Cum_time[j];
                            modec_out_.width(8);
                            modec_out_ << "s";
                        }
                    }

                    modec_out_ << '\n';
                    modec_out_ << '\n';

                    for (unsigned int i = 0; i < Nucl_size; ++i) {
                        modec_out_.width(8);
                        int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                        modec_out_ << nucl_id;
                        int NZ = nucl_id / 10000;
                        int NA = (nucl_id - NZ * 10000) / 10;
                        stringstream ss;
                        ss << NA;
                        int NG = nucl_id - nucl_id / 10 * 10;
                        string Nucl_Name;
                        Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                        if (NG == 1) {
                            Nucl_Name += "m";
                        }
                        modec_out_.width(8);
                        modec_out_ << Nucl_Name;
                        if (print_mode_ == 0) {
                            modec_out_.width(18);
                            modec_out_.setf(ios::scientific | ios::uppercase);
                            modec_out_.precision(9);
                            modec_out_ << n_vector_[0][i + Nucl_size]*Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
                            modec_out_.width(18);
                            if (n_vector_[size_tot][i + Nucl_size] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[size_tot][i + Nucl_size]*Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
                            }

                            modec_out_ << '\n';
                        }
                        if (print_mode_ == 1) {
                            modec_out_.setf(ios::scientific | ios::uppercase);
                            modec_out_.precision(9);
                            for (int j = 0; j <= size_tot; ++j) {
                                modec_out_.width(18);
                                if (n_vector_[j][i + Nucl_size] < cutoff) {
                                    modec_out_ << 0.0;
                                } else {
                                    modec_out_ << n_vector_[j][i + Nucl_size]*Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
                                }
                            }
                            modec_out_ << '\n';
                        }
                    }

                    modec_out_.close();

                }

                if (if_print_stockage_toxicity_ == true) {
                    vector<double > lamda;
                    lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

                    vector<double > TOXCI;
                    TOXCI = ModecNuclideLibrary.nuclide_library_vector_[10];

                    output_filename_ = input_filename_ + ".stockage.toxicity";
                    output_file_ = work_direc_ + output_filename_;
                    modec_out_.open(output_file_);
                    modec_out_ << "*********************************** Nuclide Radiotoxicity, Sv ***********************************";
                    modec_out_ << '\n';
                    modec_out_ << '\n';

                    modec_out_.setf(ios::left);
                    int size = evolution_mode_.size();
                    int size_tot = 0;
                    for (int i = 0; i < size; ++i) {
                        size_tot += substep_[i];
                    }

                    vector<double> Cum_time;
                    Cum_time.resize(size_tot + 1);

                    int count = 0;
                    for (int i = 0; i < size; ++i) {
                        for (int j = 0; j < substep_[i]; ++j) {
                            Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                            count++;
                        }
                    }


                    modec_out_.width(8);
                    modec_out_ << "";
                    modec_out_.width(8);
                    modec_out_ << "";
                    modec_out_.width(18);
                    modec_out_ << "Initial";

                    if (print_mode_ == 0) {
                        modec_out_.width(10);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                        modec_out_ << Cum_time[size_tot];
                        modec_out_.width(8);
                        modec_out_ << "s";
                    } else {
                        for (int j = 1; j <= size_tot; ++j) {
                            modec_out_.width(10);
                            modec_out_.setf(ios::scientific | ios::uppercase);
                            modec_out_.precision(4);
                            modec_out_ << Cum_time[j];
                            modec_out_.width(8);
                            modec_out_ << "s";
                        }
                    }

                    modec_out_ << '\n';
                    modec_out_ << '\n';

                    for (unsigned int i = 0; i < Nucl_size; ++i) {
                        modec_out_.width(8);
                        int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                        modec_out_ << nucl_id;
                        int NZ = nucl_id / 10000;
                        int NA = (nucl_id - NZ * 10000) / 10;
                        stringstream ss;
                        ss << NA;
                        int NG = nucl_id - nucl_id / 10 * 10;
                        string Nucl_Name;
                        Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                        if (NG == 1) {
                            Nucl_Name += "m";
                        }
                        modec_out_.width(8);
                        modec_out_ << Nucl_Name;
                        if (print_mode_ == 0) {
                            modec_out_.width(18);
                            modec_out_.setf(ios::scientific | ios::uppercase);
                            modec_out_.precision(9);
                            modec_out_ << n_vector_[0][i + Nucl_size]*Avogadro_Constant*lamda[i] * TOXCI[i];
                            modec_out_.width(18);
                            if (n_vector_[size_tot][i + Nucl_size] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[size_tot][i + Nucl_size]*Avogadro_Constant*lamda[i] * TOXCI[i];
                            }

                            modec_out_ << '\n';
                        }
                        if (print_mode_ == 1) {
                            modec_out_.setf(ios::scientific | ios::uppercase);
                            modec_out_.precision(9);
                            for (int j = 0; j <= size_tot; ++j) {
                                modec_out_.width(18);
                                if (n_vector_[j][i + Nucl_size] < cutoff) {
                                    modec_out_ << 0.0;
                                } else {
                                    modec_out_ << n_vector_[j][i + Nucl_size]*Avogadro_Constant*lamda[i] * TOXCI[i];
                                }
                            }
                            modec_out_ << '\n';
                        }
                    }

                    modec_out_.close();

                }
            }
        } else {
            if (input_filename_ != "modec.inp") {
                output_filename_ = input_filename_ + ".eql.concentration";
            }
            output_file_ = work_direc_ + output_filename_;
            modec_out_.open(output_file_);
            modec_out_ << "*********************************** Nuclide Equilibrium Concentrations, " << dens_unit_ << " ***********************************";
            modec_out_ << '\n';
            modec_out_ << '\n';
            modec_out_.setf(ios::left);

            modec_out_.width(8);
            modec_out_ << "";
            modec_out_.width(8);
            modec_out_ << "";
            modec_out_.width(18);
            modec_out_ << "Initial";

            modec_out_.width(18);
            modec_out_ << "Equilibrium";

            modec_out_ << '\n';
            ///////////////////////  flux_ ///////////////////
            modec_out_.width(8);
            modec_out_ << "";
            modec_out_.width(8);
            modec_out_ << "flux_: ";

            modec_out_.width(10);
            modec_out_.setf(ios::scientific | ios::uppercase);
            modec_out_.precision(4);
            modec_out_ << flux_vector_[0];
            modec_out_.width(8);
            modec_out_ << "";

            modec_out_.width(10);
            modec_out_.setf(ios::scientific | ios::uppercase);
            modec_out_.precision(4);
            modec_out_ << flux_vector_[1];
            modec_out_.width(8);
            modec_out_ << "";

            modec_out_ << '\n';

            ///////////////////////  power ///////////////////
            modec_out_.width(8);
            modec_out_ << "";
            modec_out_.width(8);
            modec_out_ << "Power: ";

            modec_out_.width(10);
            modec_out_.setf(ios::scientific | ios::uppercase);
            modec_out_.precision(4);
            modec_out_ << power_vector_[0];
            modec_out_.width(8);
            modec_out_ << "";

            modec_out_.width(10);
            modec_out_.setf(ios::scientific | ios::uppercase);
            modec_out_.precision(4);
            modec_out_ << power_vector_[1];
            modec_out_.width(8);
            modec_out_ << "";
            modec_out_ << '\n';

            modec_out_ << '\n';
            ////////////////////////////////////////////////////////

            for (unsigned int i = 0; i < Nucl_size; ++i) {
                modec_out_.width(8);
                int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                modec_out_ << nucl_id;
                int NZ = nucl_id / 10000;
                int NA = (nucl_id - NZ * 10000) / 10;
                stringstream ss;
                ss << NA;
                int NG = nucl_id - nucl_id / 10 * 10;
                string Nucl_Name;
                Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                if (NG == 1) {
                    Nucl_Name += "m";
                }
                modec_out_.width(8);
                modec_out_ << Nucl_Name;

                double convert_coeff;

                if (dens_unit_ == "mol") {
                    convert_coeff = 1.0;
                } else if (dens_unit_ == "g") {
                    convert_coeff = double(NA);
                } else if (dens_unit_ == "kg") {
                    convert_coeff = double(NA) / 1000.0;
                } else if (dens_unit_ == "atom") {
                    convert_coeff = Avogadro_Constant;
                } else if (dens_unit_ == "atom/(barn-cm)") {
                    convert_coeff = Avogadro_Constant * (1.0e-24);
                }

                modec_out_.width(18);
                modec_out_.setf(ios::scientific | ios::uppercase);
                modec_out_.precision(9);
                modec_out_ << n_vector_[0][i] * convert_coeff;
                modec_out_.width(18);
                if (n_vector_[1][i] < cutoff) {
                    modec_out_ << 0.0;
                } else {
                    modec_out_ << n_vector_[1][i] * convert_coeff;
                }

                modec_out_ << '\n';

            }

            modec_out_.close();

            if (if_print_activity_ == true) {

                vector<double > lamda;
                lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

                output_filename_ = input_filename_ + ".eql.radioactivity";
                output_file_ = work_direc_ + output_filename_;
                modec_out_.open(output_file_);
                modec_out_ << "*********************************** Nuclide Equilibrium Radioactivity, Ci ***********************************";
                modec_out_ << '\n';
                modec_out_ << '\n';
                modec_out_.setf(ios::left);

                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(18);
                modec_out_ << "Initial";

                modec_out_.width(18);
                modec_out_ << "Equilibrium";

                modec_out_ << '\n';

                modec_out_ << '\n';
                modec_out_ << '\n';

                for (unsigned int i = 0; i < Nucl_size; ++i) {
                    modec_out_.width(8);
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    modec_out_ << nucl_id;
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    stringstream ss;
                    ss << NA;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string Nucl_Name;
                    Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                    if (NG == 1) {
                        Nucl_Name += "m";
                    }
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;

                    modec_out_.width(18);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(9);
                    modec_out_ << n_vector_[0][i]*Avogadro_Constant*lamda[i] / (3.7e10);
                    modec_out_.width(18);
                    if (n_vector_[1][i] < cutoff) {
                        modec_out_ << 0.0;
                    } else {
                        modec_out_ << n_vector_[1][i]*Avogadro_Constant*lamda[i] / (3.7e10); // 单位为Ci
                    }

                    modec_out_ << '\n';

                }

                modec_out_.close();

            }

            if (if_print_decayenergy_ == true) {
                vector<double > lamda;
                lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

                vector<double > Q;
                Q = ModecNuclideLibrary.nuclide_library_vector_[3];

                output_filename_ = input_filename_ + ".eql.decayheat";
                output_file_ = work_direc_ + output_filename_;
                modec_out_.open(output_file_);
                modec_out_ << "*********************************** Nuclide Equilibrium Decay Heat, Watts ***********************************";
                modec_out_ << '\n';
                modec_out_ << '\n';
                modec_out_.setf(ios::left);

                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(18);
                modec_out_ << "Initial";

                modec_out_.width(18);
                modec_out_ << "Equilibrium";

                modec_out_ << '\n';

                modec_out_ << '\n';
                modec_out_ << '\n';


                for (unsigned int i = 0; i < Nucl_size; ++i) {
                    modec_out_.width(8);
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    modec_out_ << nucl_id;
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    stringstream ss;
                    ss << NA;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string Nucl_Name;
                    Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                    if (NG == 1) {
                        Nucl_Name += "m";
                    }
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;

                    modec_out_.width(18);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(9);
                    modec_out_ << n_vector_[0][i]*Avogadro_Constant*lamda[i] * Q[i] * Electron_Coulomb*1.0e6;
                    modec_out_.width(18);
                    if (n_vector_[1][i] < cutoff) {
                        modec_out_ << 0.0;
                    } else {
                        modec_out_ << n_vector_[1][i]*Avogadro_Constant*lamda[i] * Q[i] * Electron_Coulomb*1.0e6; // 转换成单位 W
                    }

                    modec_out_ << '\n';
                }

                modec_out_.close();

            }

            if (if_print_ampc_ == true) {
                vector<double > lamda;
                lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

                vector<double > AMPC;
                AMPC = ModecNuclideLibrary.nuclide_library_vector_[4];

                output_filename_ = input_filename_ + ".eql.ampc";
                output_file_ = work_direc_ + output_filename_;
                modec_out_.open(output_file_);
                modec_out_ << "*********************************** Nuclide Equilibrium Radioactivity, m3-air/Bq ***********************************";
                modec_out_ << '\n';
                modec_out_ << '\n';
                modec_out_.setf(ios::left);

                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(18);
                modec_out_ << "Initial";

                modec_out_.width(18);
                modec_out_ << "Equilibrium";

                modec_out_ << '\n';

                modec_out_ << '\n';
                modec_out_ << '\n';


                for (unsigned int i = 0; i < Nucl_size; ++i) {
                    modec_out_.width(8);
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    modec_out_ << nucl_id;
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    stringstream ss;
                    ss << NA;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string Nucl_Name;
                    Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                    if (NG == 1) {
                        Nucl_Name += "m";
                    }
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;

                    modec_out_.width(18);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(9);
                    modec_out_ << n_vector_[0][i]*Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
                    modec_out_.width(18);
                    if (n_vector_[1][i] < cutoff) {
                        modec_out_ << 0.0;
                    } else {
                        modec_out_ << n_vector_[1][i]*Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
                    }

                    modec_out_ << '\n';

                }

                modec_out_.close();

            }

            if (if_print_wmpc_ == true) {
                vector<double > lamda;
                lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

                vector<double > WMPC;
                WMPC = ModecNuclideLibrary.nuclide_library_vector_[4];

                output_filename_ = input_filename_ + ".eql.wmpc";
                output_file_ = work_direc_ + output_filename_;
                modec_out_.open(output_file_);
                modec_out_ << "*********************************** Nuclide Equilibrium Radioactivity, m3-air/Bq ***********************************";
                modec_out_ << '\n';
                modec_out_ << '\n';
                modec_out_.setf(ios::left);

                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(18);
                modec_out_ << "Initial";

                modec_out_.width(18);
                modec_out_ << "Equilibrium";

                modec_out_ << '\n';

                modec_out_ << '\n';
                modec_out_ << '\n';


                for (unsigned int i = 0; i < Nucl_size; ++i) {
                    modec_out_.width(8);
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    modec_out_ << nucl_id;
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    stringstream ss;
                    ss << NA;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string Nucl_Name;
                    Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                    if (NG == 1) {
                        Nucl_Name += "m";
                    }
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;

                    modec_out_.width(18);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(9);
                    modec_out_ << n_vector_[0][i]*Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
                    modec_out_.width(18);
                    if (n_vector_[1][i] < cutoff) {
                        modec_out_ << 0.0;
                    } else {
                        modec_out_ << n_vector_[1][i]*Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
                    }

                    modec_out_ << '\n';

                }

                modec_out_.close();

            }

            if (if_print_toxicity_ == true) {
                vector<double > lamda;
                lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

                vector<double > TOXIC;
                TOXIC = ModecNuclideLibrary.nuclide_library_vector_[10];

                output_filename_ = input_filename_ + ".eql.toxicity";
                output_file_ = work_direc_ + output_filename_;
                modec_out_.open(output_file_);
                modec_out_ << "*********************************** Nuclide Equilibrium Radiotoxicity, Sv ***********************************";
                modec_out_ << '\n';
                modec_out_ << '\n';
                modec_out_.setf(ios::left);

                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(18);
                modec_out_ << "Initial";

                modec_out_.width(18);
                modec_out_ << "Equilibrium";

                modec_out_ << '\n';

                modec_out_ << '\n';
                modec_out_ << '\n';


                for (unsigned int i = 0; i < Nucl_size; ++i) {
                    modec_out_.width(8);
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    modec_out_ << nucl_id;
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    stringstream ss;
                    ss << NA;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string Nucl_Name;
                    Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                    if (NG == 1) {
                        Nucl_Name += "m";
                    }
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;

                    modec_out_.width(18);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(9);
                    modec_out_ << n_vector_[0][i]*Avogadro_Constant*lamda[i] * TOXIC[i];
                    modec_out_.width(18);
                    if (n_vector_[1][i] < cutoff) {
                        modec_out_ << 0.0;
                    } else {
                        modec_out_ << n_vector_[1][i]*Avogadro_Constant*lamda[i] * TOXIC[i];
                    }

                    modec_out_ << '\n';

                }

                modec_out_.close();

            }
        }
    } else { // 流动情况下的结果输出
        if (input_filename_ != "modec.inp") {
            output_filename_ = input_filename_ + ".concentration";
        }

        output_file_ = work_direc_ + output_filename_;
        modec_out_.open(output_file_);
        modec_out_ << "*********************************** Nuclide Concentrations, " << dens_unit_ << " ***********************************";
        modec_out_ << '\n';
        modec_out_ << '\n';
        modec_out_.setf(ios::left);
        int size = evolution_mode_.size();
        int size_tot = 0;
        vector<string> time_unit;
        for (int i = 0; i < size; ++i) {
            size_tot += substep_[i];
        }

        vector<double> Cum_time;
        Cum_time.resize(size_tot + 1);

        int count = 0;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < substep_[i]; ++j) {
                Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                time_unit.push_back(time_unit_[i]);
                count++;
            }
        }

        if (print_mode_ == 0) {
            modec_out_.width(10);
            modec_out_.setf(ios::right);
            if (time_unit[size_tot - 1] == "d" || time_unit[size_tot - 1] == "day" || time_unit[size_tot - 1] == "days") {
                Cum_time[size_tot] = Cum_time[size_tot] / 24.0 / 3600.0;
            }
            if (time_unit[size_tot - 1] == "m" || time_unit[size_tot - 1] == "month" || time_unit[size_tot - 1] == "months") {
                Cum_time[size_tot] = Cum_time[size_tot] / 24.0 / 3600.0 / 30.0;
            } else if (time_unit[size_tot - 1] == "y" || time_unit[size_tot - 1] == "year" || time_unit[size_tot - 1] == "years") {
                Cum_time[size_tot] = Cum_time[size_tot] / 24.0 / 3600.0 / 365.25;
                time_unit[size_tot - 1] = "yr";
            }

            if (Cum_time[size_tot] >= 1.0e+3) {
                modec_out_.setf(ios::scientific | ios::uppercase);
                modec_out_.precision(4);
            } else {
                modec_out_.precision(5);
            }

            modec_out_ << Cum_time[size_tot];
            modec_out_.unsetf(ios::right);
            modec_out_ << " ";
            modec_out_.width(7);
            modec_out_ << time_unit[size_tot - 1];
        } else {
            for (int j = 1; j <= size_tot; ++j) {
                modec_out_.setf(ios::right);
                if (time_unit[j - 1] == "d" || time_unit[j - 1] == "day" || time_unit[j - 1] == "days") {
                    Cum_time[j] = Cum_time[j] / 24.0 / 3600.0;
                }
                if (time_unit[j - 1] == "m" || time_unit[j - 1] == "month" || time_unit[j - 1] == "months") {
                    Cum_time[j] = Cum_time[j] / 24.0 / 3600.0 / 30.0;
                } else if (time_unit[j - 1] == "y" || time_unit[j - 1] == "year" || time_unit[j - 1] == "years") {
                    Cum_time[j] = Cum_time[j] / 24.0 / 3600.0 / 365.25;
                    time_unit[j - 1] = "yr";
                }
                modec_out_.width(10);
                if (Cum_time[j] >= 1.0e+3) {
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                } else {
                    modec_out_.precision(5);
                }
                modec_out_ << Cum_time[j];
                modec_out_.unsetf(ios::right);
                modec_out_ << " ";
                modec_out_.width(7);
                modec_out_ << time_unit[j - 1];
            }
        }
        modec_out_ << '\n';

        modec_out_.width(8);
        modec_out_ << "";
        modec_out_.width(8);
        modec_out_ << "Power: ";

        if (print_mode_ == 0) {
            modec_out_.width(10);
            modec_out_.setf(ios::scientific | ios::uppercase);
            modec_out_.precision(4);
            modec_out_ << power_vector_[0];
            modec_out_.width(8);
            modec_out_ << "MW";

            modec_out_.width(10);
            modec_out_.setf(ios::scientific | ios::uppercase);
            modec_out_.precision(4);
            modec_out_ << power_vector_[size_tot];
            modec_out_.width(8);
            modec_out_ << "MW";
        } else {
            /*modec_out_.width(10);
            modec_out_.setf(ios::scientific | ios::uppercase);
            modec_out_.precision(4);
            modec_out_ << 0.0;
            modec_out_.width(8);
            modec_out_ << "MW";*/

            for (int j = 0; j <= size_tot; ++j) {
                modec_out_.width(10);
                modec_out_.setf(ios::scientific | ios::uppercase);
                modec_out_.precision(4);
                modec_out_ << power_vector_[j];
                modec_out_.width(8);
                modec_out_ << " MW";
            }
        }

        modec_out_ << '\n';

        modec_out_.width(8);
        modec_out_ << "";
        modec_out_.width(8);
        modec_out_ << "Flux: ";

        if (print_mode_ == 0) {
            modec_out_.width(10);
            modec_out_.setf(ios::scientific | ios::uppercase);
            modec_out_.precision(4);
            modec_out_ << flux_vector_[0];
            modec_out_.width(8);
            modec_out_ << "";

            modec_out_.width(10);
            modec_out_.setf(ios::scientific | ios::uppercase);
            modec_out_.precision(4);
            modec_out_ << flux_vector_[size_tot];
            modec_out_.width(8);
            modec_out_ << "";
        } else {
            //modec_out_.width(10);
            //modec_out_.setf(ios::scientific | ios::uppercase);
            //modec_out_.precision(4);
            //modec_out_ << 0.0;
            //modec_out_.width(8);
            //modec_out_ << "";
            for (int j = 0; j <= size_tot; ++j) {
                modec_out_.width(10);
                modec_out_.setf(ios::scientific | ios::uppercase);
                modec_out_.precision(4);
                modec_out_ << flux_vector_[j];
                modec_out_.width(8);
                modec_out_ << "";
            }
        }
        modec_out_ << '\n';

        //////////////    Calculate kinf_vector_, Neutron_Production_Rate and Neutron_Absorption_Rate      /////////
        vector<double> kinf;
        vector<double> prod_neu;
        vector<double> absorption_neu;
        if (if_print_kinf_ != 0 || if_print_fission_rate_ != 0 || if_print_absorption_rate_ != 0) {
            kinf.resize(size_tot + 1);
            prod_neu.resize(size_tot + 1);
            absorption_neu.resize(size_tot + 1);

            for (int j = 0; j <= size_tot; ++j) {
                for (int i = 0; i < Nucl_size; ++i) {
                    prod_neu[j] += n_vector_[j][2 * i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24;
                    absorption_neu[j] += n_vector_[j][2 * i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24;
                }
                prod_neu[j] = prod_neu[j] * flux_vector_[j];
                absorption_neu[j] = absorption_neu[j] * flux_vector_[j];

                if (absorption_neu[j] != 0.0) {
                    kinf[j] = prod_neu[j] / absorption_neu[j];
                }

            }

            if (if_print_kinf_ == true) {

                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "kinf: ";

                if (print_mode_ == 0) {
                    modec_out_.width(18);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(9);
                    modec_out_ << kinf[0];
                    modec_out_.width(18);
                    modec_out_ << kinf[size_tot];
                    modec_out_ << '\n';
                }
                if (print_mode_ == 1) {
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(9);
                    for (int j = 0; j <= size_tot; ++j) {
                        modec_out_.width(18);
                        modec_out_ << kinf[j];
                    }
                    modec_out_ << '\n';
                }
            }

            if (if_print_fission_rate_ != 0) {

                modec_out_.width(2);
                modec_out_ << "";
                modec_out_.width(14);
                modec_out_ << "NeuProdRate: ";

                if (print_mode_ == 0) {
                    modec_out_.width(18);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(9);
                    modec_out_ << prod_neu[0];
                    modec_out_.width(18);
                    modec_out_ << prod_neu[size_tot];
                    modec_out_ << '\n';
                }
                if (print_mode_ == 1) {
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(9);
                    for (int j = 0; j <= size_tot; ++j) {
                        modec_out_.width(18);
                        modec_out_ << prod_neu[j];
                    }
                    modec_out_ << '\n';
                }
            }

            if (if_print_absorption_rate_ != 0) {

                modec_out_.width(2);
                modec_out_ << "";
                modec_out_.width(14);
                modec_out_ << "NeuAbsRate: ";

                if (print_mode_ == 0) {
                    modec_out_.width(18);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(9);
                    modec_out_ << absorption_neu[0];
                    modec_out_.width(18);
                    modec_out_ << absorption_neu[size_tot];
                    modec_out_ << '\n';
                }
                if (print_mode_ == 1) {
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(9);
                    for (int j = 0; j <= size_tot; ++j) {
                        modec_out_.width(18);
                        modec_out_ << absorption_neu[j];
                    }
                    modec_out_ << '\n';
                }
            }
            modec_out_ << '\n';
        }
        //////////////////////////////////////////////////////////////////////////////////////////////////

        modec_out_ << '\n';

        for (unsigned int i = 0; i < Nucl_size; ++i) {
            modec_out_.width(8);
            int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
            modec_out_ << nucl_id;
            int NZ = nucl_id / 10000;
            int NA = (nucl_id - NZ * 10000) / 10;
            stringstream ss;
            ss << NA;
            int NG = nucl_id - nucl_id / 10 * 10;
            string Nucl_Name;
            Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
            if (NG == 1) {
                Nucl_Name += "m";
            }
            modec_out_.width(8);
            modec_out_ << Nucl_Name;

            double convert_coeff;

            if (dens_unit_ == "mol") {
                convert_coeff = 1.0;
            } else if (dens_unit_ == "g") {
                convert_coeff = double(NA);
            } else if (dens_unit_ == "kg") {
                convert_coeff = double(NA) / 1000.0;
            } else if (dens_unit_ == "atom") {
                convert_coeff = Avogadro_Constant;
            } else if (dens_unit_ == "atom/(barn-cm)") {
                convert_coeff = Avogadro_Constant * (1.0e-24);
            }

            if (print_mode_ == 0) {
                modec_out_.width(18);
                modec_out_.setf(ios::scientific | ios::uppercase);
                modec_out_.precision(9);
                modec_out_ << n_vector_[0][2 * i] * convert_coeff;
                modec_out_.width(18);
                if (n_vector_[size_tot][2 * i] < cutoff) {
                    modec_out_ << 0.0;
                } else {
                    modec_out_ << n_vector_[size_tot][2 * i] * convert_coeff;
                }

                modec_out_ << '\n';

                modec_out_.width(8);
                modec_out_ << nucl_id;
                modec_out_.width(8);
                modec_out_ << Nucl_Name;
                modec_out_.width(18);
                modec_out_.setf(ios::scientific | ios::uppercase);
                modec_out_.precision(9);
                modec_out_ << n_vector_[0][2 * i + 1] * convert_coeff;
                modec_out_.width(18);
                if (n_vector_[size_tot][2 * i + 1] < cutoff) {
                    modec_out_ << 0.0;
                } else {
                    modec_out_ << n_vector_[size_tot][2 * i + 1] * convert_coeff;
                }

                modec_out_ << '\n';
            }
            if (print_mode_ == 1) {
                modec_out_.setf(ios::scientific | ios::uppercase);
                modec_out_.precision(9);
                for (int j = 0; j <= size_tot; ++j) {
                    modec_out_.width(18);
                    if (n_vector_[j][2 * i] < cutoff ) {
                        modec_out_ << 0.0;
                    } else {
                        modec_out_ << n_vector_[j][2 * i] * convert_coeff;
                    }
                }
                modec_out_ << '\n';

                modec_out_.width(8);
                modec_out_ << nucl_id;
                modec_out_.width(8);
                modec_out_ << Nucl_Name;
                modec_out_.setf(ios::scientific | ios::uppercase);
                modec_out_.precision(9);
                for (int j = 0; j <= size_tot; ++j) {
                    modec_out_.width(18);
                    if (n_vector_[j][2 * i + 1] < cutoff) {
                        modec_out_ << 0.0;
                    } else {
                        modec_out_ << n_vector_[j][2 * i + 1] * convert_coeff;
                    }
                }
                modec_out_ << '\n';
            }
        }

        modec_out_.close();

        if (if_print_activity_ == true) {

            vector<double > lamda;
            lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

            output_filename_ = input_filename_ + ".radioactivity";

            output_file_ = work_direc_ + output_filename_;
            modec_out_.open(output_file_);
            modec_out_ << "*********************************** Nuclide Radioactivity, Ci ***********************************";
            modec_out_ << '\n';
            modec_out_ << '\n';

            modec_out_.setf(ios::left);
            int size = evolution_mode_.size();
            int size_tot = 0;
            for (int i = 0; i < size; ++i) {
                size_tot += substep_[i];
            }

            vector<double> Cum_time;
            Cum_time.resize(size_tot + 1);

            int count = 0;
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < substep_[i]; ++j) {
                    Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                    count++;
                }
            }


            modec_out_.width(8);
            modec_out_ << "";
            modec_out_.width(8);
            modec_out_ << "";
            modec_out_.width(18);
            modec_out_ << "Initial";

            if (print_mode_ == 0) {
                modec_out_.width(10);
                modec_out_.setf(ios::scientific | ios::uppercase);
                modec_out_.precision(4);
                modec_out_ << Cum_time[size_tot];
                modec_out_.width(8);
                modec_out_ << "s";
            } else {
                for (int j = 1; j <= size_tot; ++j) {
                    modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << Cum_time[j];
                    modec_out_.width(8);
                    modec_out_ << "s";
                }
            }

            modec_out_ << '\n';
            modec_out_ << '\n';

            for (unsigned int i = 0; i < Nucl_size; ++i) {
                modec_out_.width(8);
                int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                modec_out_ << nucl_id;
                int NZ = nucl_id / 10000;
                int NA = (nucl_id - NZ * 10000) / 10;
                stringstream ss;
                ss << NA;
                int NG = nucl_id - nucl_id / 10 * 10;
                string Nucl_Name;
                Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                if (NG == 1) {
                    Nucl_Name += "m";
                }
                modec_out_.width(8);
                modec_out_ << Nucl_Name;
                if (print_mode_ == 0) {
                    modec_out_.width(18);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(9);
                    modec_out_ << n_vector_[0][2 * i] * Avogadro_Constant*lamda[i] / (3.7e10);
                    modec_out_.width(18);
                    if (n_vector_[size_tot][2 * i] < cutoff) {
                        modec_out_ << 0.0;
                    } else {
                        modec_out_ << n_vector_[size_tot][2 * i] * Avogadro_Constant*lamda[i] / (3.7e10); // 单位为Ci
                    }

                    modec_out_ << '\n';

                    modec_out_.width(8);
                    modec_out_ << nucl_id;
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;
                    modec_out_.width(18);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(9);
                    modec_out_ << n_vector_[0][2 * i + 1] * Avogadro_Constant*lamda[i] / (3.7e10);
                    modec_out_.width(18);
                    if (n_vector_[size_tot][2 * i + 1] < cutoff) {
                        modec_out_ << 0.0;
                    } else {
                        modec_out_ << n_vector_[size_tot][2 * i + 1] * Avogadro_Constant*lamda[i] / (3.7e10); // 单位为Ci
                    }

                    modec_out_ << '\n';

                }
                if (print_mode_ == 1) {
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(9);
                    for (int j = 0; j <= size_tot; ++j) {
                        modec_out_.width(18);
                        if (n_vector_[j][2 * i] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[j][2 * i] * Avogadro_Constant*lamda[i] / (3.7e10); // 单位为Ci
                        }
                    }
                    modec_out_ << '\n';

                    modec_out_.width(8);
                    modec_out_ << nucl_id;
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(9);
                    for (int j = 0; j <= size_tot; ++j) {
                        modec_out_.width(18);
                        if (n_vector_[j][2 * i + 1] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[j][2 * i + 1] * Avogadro_Constant*lamda[i] / (3.7e10); // 单位为Ci
                        }
                    }
                    modec_out_ << '\n';
                }
            }

            modec_out_.close();

        }

        if (if_print_decayenergy_ == true) {
            vector<double > lamda;
            lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

            vector<double > Q;
            Q = ModecNuclideLibrary.nuclide_library_vector_[3];

            output_filename_ = input_filename_ + ".decayheat";

            output_file_ = work_direc_ + output_filename_;
            modec_out_.open(output_file_);
            modec_out_ << "*********************************** Nuclide Decay Heat, Watts ***********************************";
            modec_out_ << '\n';
            modec_out_ << '\n';

            modec_out_.setf(ios::left);
            int size = evolution_mode_.size();
            int size_tot = 0;
            for (int i = 0; i < size; ++i) {
                size_tot += substep_[i];
            }

            vector<double> Cum_time;
            Cum_time.resize(size_tot + 1);

            int count = 0;
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < substep_[i]; ++j) {
                    Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                    count++;
                }
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(18);
                modec_out_ << "Initial";

                if (print_mode_ == 0) {
                    modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << Cum_time[size_tot];
                    modec_out_.width(8);
                    modec_out_ << "s";
                } else {
                    for (int j = 1; j <= size_tot; ++j) {
                        modec_out_.width(10);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                        modec_out_ << Cum_time[j];
                        modec_out_.width(8);
                        modec_out_ << "s";
                    }
                }

                modec_out_ << '\n';
                modec_out_ << '\n';

                for (unsigned int i = 0; i < Nucl_size; ++i) {
                    modec_out_.width(8);
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    modec_out_ << nucl_id;
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    stringstream ss;
                    ss << NA;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string Nucl_Name;
                    Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                    if (NG == 1) {
                        Nucl_Name += "m";
                    }
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;
                    if (print_mode_ == 0) {
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << n_vector_[0][2 * i] * Avogadro_Constant*lamda[i] * Q[i] * Electron_Coulomb*1.0e6;
                        modec_out_.width(18);
                        if (n_vector_[size_tot][2 * i] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[size_tot][2 * i] * Avogadro_Constant*lamda[i] * Q[i] * Electron_Coulomb*1.0e6; // 转换成单位 W
                        }

                        modec_out_ << '\n';

                        modec_out_.width(8);
                        modec_out_ << nucl_id;
                        modec_out_.width(8);
                        modec_out_ << Nucl_Name;
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << n_vector_[0][2 * i + 1] * Avogadro_Constant*lamda[i] * Q[i] * Electron_Coulomb*1.0e6;
                        modec_out_.width(18);
                        if (n_vector_[size_tot][2 * i + 1] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[size_tot][2 * i + 1] * Avogadro_Constant*lamda[i] * Q[i] * Electron_Coulomb*1.0e6; // 转换成单位 W
                        }

                        modec_out_ << '\n';

                    }
                    if (print_mode_ == 1) {
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            if (n_vector_[j][2 * i] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[j][2 * i] * Avogadro_Constant*lamda[i] * Q[i] * Electron_Coulomb*1.0e6;
                            }
                        }
                        modec_out_ << '\n';

                        modec_out_.width(8);
                        modec_out_ << nucl_id;
                        modec_out_.width(8);
                        modec_out_ << Nucl_Name;
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            if (n_vector_[j][2 * i + 1] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[j][2 * i + 1] * Avogadro_Constant*lamda[i] * Q[i] * Electron_Coulomb*1.0e6;
                            }
                        }
                        modec_out_ << '\n';
                    }
                }

                modec_out_.close();

            }

            if (if_print_ampc_ == true) {
                vector<double > lamda;
                lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

                vector<double > AMPC;
                AMPC = ModecNuclideLibrary.nuclide_library_vector_[4];

                output_filename_ = input_filename_ + ".ampc";

                output_file_ = work_direc_ + output_filename_;
                modec_out_.open(output_file_);
                modec_out_ << "*********************************** Nuclide Radioactivity, m3-air/Bq ***********************************";
                modec_out_ << '\n';
                modec_out_ << '\n';

                modec_out_.setf(ios::left);
                int size = evolution_mode_.size();
                int size_tot = 0;
                for (int i = 0; i < size; ++i) {
                    size_tot += substep_[i];
                }

                vector<double> Cum_time;
                Cum_time.resize(size_tot + 1);

                int count = 0;
                for (int i = 0; i < size; ++i) {
                    for (int j = 0; j < substep_[i]; ++j) {
                        Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                        count++;
                    }
                }


                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(18);
                modec_out_ << "Initial";

                if (print_mode_ == 0) {
                    modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << Cum_time[size_tot];
                    modec_out_.width(8);
                    modec_out_ << "s";
                } else {
                    for (int j = 1; j <= size_tot; ++j) {
                        modec_out_.width(10);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                        modec_out_ << Cum_time[j];
                        modec_out_.width(8);
                        modec_out_ << "s";
                    }
                }

                modec_out_ << '\n';
                modec_out_ << '\n';

                for (unsigned int i = 0; i < Nucl_size; ++i) {
                    modec_out_.width(8);
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    modec_out_ << nucl_id;
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    stringstream ss;
                    ss << NA;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string Nucl_Name;
                    Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                    if (NG == 1) {
                        Nucl_Name += "m";
                    }
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;
                    if (print_mode_ == 0) {
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << n_vector_[0][2 * i] * Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
                        modec_out_.width(18);
                        if (n_vector_[size_tot][2 * i] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[size_tot][2 * i] * Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
                        }

                        modec_out_ << '\n';

                        modec_out_.width(8);
                        modec_out_ << nucl_id;
                        modec_out_.width(8);
                        modec_out_ << Nucl_Name;
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << n_vector_[0][2 * i + 1] * Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
                        modec_out_.width(18);
                        if (n_vector_[size_tot][2 * i + 1] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[size_tot][2 * i + 1] * Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
                        }

                        modec_out_ << '\n';
                    }
                    if (print_mode_ == 1) {
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            if (n_vector_[j][2 * i] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[j][2 * i] * Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
                            }
                        }
                        modec_out_ << '\n';

                        modec_out_.width(8);
                        modec_out_ << nucl_id;
                        modec_out_.width(8);
                        modec_out_ << Nucl_Name;
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            if (n_vector_[j][2 * i + 1] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[j][2 * i + 1] * Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
                            }
                        }
                        modec_out_ << '\n';
                    }
                }

                modec_out_.close();

            }

            if (if_print_wmpc_ == true) {
                vector<double > lamda;
                lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

                vector<double > WMPC;
                WMPC = ModecNuclideLibrary.nuclide_library_vector_[4];

                output_filename_ = input_filename_ + ".wmpc";

                output_file_ = work_direc_ + output_filename_;
                modec_out_.open(output_file_);
                modec_out_ << "*********************************** Nuclide Radioactivity, m3-air/Bq ***********************************";
                modec_out_ << '\n';
                modec_out_ << '\n';

                modec_out_.setf(ios::left);
                int size = evolution_mode_.size();
                int size_tot = 0;
                for (int i = 0; i < size; ++i) {
                    size_tot += substep_[i];
                }

                vector<double> Cum_time;
                Cum_time.resize(size_tot + 1);

                int count = 0;
                for (int i = 0; i < size; ++i) {
                    for (int j = 0; j < substep_[i]; ++j) {
                        Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                        count++;
                    }
                }


                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(18);
                modec_out_ << "Initial";

                if (print_mode_ == 0) {
                    modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << Cum_time[size_tot];
                    modec_out_.width(8);
                    modec_out_ << "s";
                } else {
                    for (int j = 1; j <= size_tot; ++j) {
                        modec_out_.width(10);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                        modec_out_ << Cum_time[j];
                        modec_out_.width(8);
                        modec_out_ << "s";
                    }
                }

                modec_out_ << '\n';
                modec_out_ << '\n';

                for (unsigned int i = 0; i < Nucl_size; ++i) {
                    modec_out_.width(8);
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    modec_out_ << nucl_id;
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    stringstream ss;
                    ss << NA;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string Nucl_Name;
                    Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                    if (NG == 1) {
                        Nucl_Name += "m";
                    }
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;
                    if (print_mode_ == 0) {
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << n_vector_[0][2 * i] * Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
                        modec_out_.width(18);
                        if (n_vector_[size_tot][2 * i] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[size_tot][2 * i] * Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
                        }

                        modec_out_ << '\n';

                        modec_out_.width(8);
                        modec_out_ << nucl_id;
                        modec_out_.width(8);
                        modec_out_ << Nucl_Name;
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << n_vector_[0][2 * i + 1] * Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
                        modec_out_.width(18);
                        if (n_vector_[size_tot][2 * i + 1] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[size_tot][2 * i + 1] * Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
                        }

                        modec_out_ << '\n';
                    }
                    if (print_mode_ == 1) {
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            if (n_vector_[j][2 * i] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[j][2 * i] * Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
                            }
                        }
                        modec_out_ << '\n';

                        modec_out_.width(8);
                        modec_out_ << nucl_id;
                        modec_out_.width(8);
                        modec_out_ << Nucl_Name;
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            if (n_vector_[j][2 * i + 1] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[j][2 * i + 1] * Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
                            }
                        }
                        modec_out_ << '\n';
                    }
                }

                modec_out_.close();

            }

            if (if_print_toxicity_ == true) {
                vector<double > lamda;
                lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

                vector<double > TOXCI;
                TOXCI = ModecNuclideLibrary.nuclide_library_vector_[10];

                output_filename_ = input_filename_ + ".toxicity";

                output_file_ = work_direc_ + output_filename_;
                modec_out_.open(output_file_);
                modec_out_ << "*********************************** Nuclide Radiotoxicity, Sv ***********************************";
                modec_out_ << '\n';
                modec_out_ << '\n';

                modec_out_.setf(ios::left);
                int size = evolution_mode_.size();
                int size_tot = 0;
                for (int i = 0; i < size; ++i) {
                    size_tot += substep_[i];
                }

                vector<double> Cum_time;
                Cum_time.resize(size_tot + 1);

                int count = 0;
                for (int i = 0; i < size; ++i) {
                    for (int j = 0; j < substep_[i]; ++j) {
                        Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                        count++;
                    }
                }


                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(18);
                modec_out_ << "Initial";

                if (print_mode_ == 0) {
                    modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << Cum_time[size_tot];
                    modec_out_.width(8);
                    modec_out_ << "s";
                } else {
                    for (int j = 1; j <= size_tot; ++j) {
                        modec_out_.width(10);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                        modec_out_ << Cum_time[j];
                        modec_out_.width(8);
                        modec_out_ << "s";
                    }
                }

                modec_out_ << '\n';
                modec_out_ << '\n';

                for (unsigned int i = 0; i < Nucl_size; ++i) {
                    modec_out_.width(8);
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    modec_out_ << nucl_id;
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    stringstream ss;
                    ss << NA;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string Nucl_Name;
                    Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                    if (NG == 1) {
                        Nucl_Name += "m";
                    }
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;
                    if (print_mode_ == 0) {
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << n_vector_[0][2 * i] * Avogadro_Constant*lamda[i] * TOXCI[i];
                        modec_out_.width(18);
                        if (n_vector_[size_tot][2 * i] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[size_tot][2 * i] * Avogadro_Constant*lamda[i] * TOXCI[i];
                        }

                        modec_out_ << '\n';

                        modec_out_.width(8);
                        modec_out_ << nucl_id;
                        modec_out_.width(8);
                        modec_out_ << Nucl_Name;
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << n_vector_[0][2 * i + 1] * Avogadro_Constant*lamda[i] * TOXCI[i];
                        modec_out_.width(18);
                        if (n_vector_[size_tot][2 * i + 1] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[size_tot][2 * i + 1] * Avogadro_Constant*lamda[i] * TOXCI[i];
                        }

                        modec_out_ << '\n';

                    }
                    if (print_mode_ == 1) {
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            if (n_vector_[j][2 * i] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[j][2 * i] * Avogadro_Constant*lamda[i] * TOXCI[i];
                            }
                        }
                        modec_out_ << '\n';

                        modec_out_.width(8);
                        modec_out_ << nucl_id;
                        modec_out_.width(8);
                        modec_out_ << Nucl_Name;
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            if (n_vector_[j][2 * i + 1] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[j][2 * i + 1] * Avogadro_Constant*lamda[i] * TOXCI[i];
                            }
                        }
                        modec_out_ << '\n';
                    }
                }

                modec_out_.close();

            }

            if (if_print_fission_rate_ == 2) {
                output_filename_ = input_filename_ + ".NeuProdRate";

                output_file_ = work_direc_ + output_filename_;
                modec_out_.open(output_file_);
                modec_out_ << "*********************************** Neutron Production Rate, s-1***********************************";
                modec_out_ << '\n';
                modec_out_ << '\n';

                modec_out_.setf(ios::left);
                int size = evolution_mode_.size();
                int size_tot = 0;
                for (int i = 0; i < size; ++i) {
                    size_tot += substep_[i];
                }

                vector<double> Cum_time;
                Cum_time.resize(size_tot + 1);

                int count = 0;
                for (int i = 0; i < size; ++i) {
                    for (int j = 0; j < substep_[i]; ++j) {
                        Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                        count++;
                    }
                }


                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(18);
                modec_out_ << "Initial";

                if (print_mode_ == 0) {
                    modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << Cum_time[size_tot];
                    modec_out_.width(8);
                    modec_out_ << "s";
                } else {
                    for (int j = 1; j <= size_tot; ++j) {
                        modec_out_.width(10);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                        modec_out_ << Cum_time[j];
                        modec_out_.width(8);
                        modec_out_ << "s";
                    }
                }

                modec_out_ << '\n';
                modec_out_ << '\n';

                for (unsigned int i = 0; i < Nucl_size; ++i) {
                    modec_out_.width(8);
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    modec_out_ << nucl_id;
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    stringstream ss;
                    ss << NA;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string Nucl_Name;
                    Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                    if (NG == 1) {
                        Nucl_Name += "m";
                    }
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;
                    if (print_mode_ == 0) {
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << n_vector_[0][2 * i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24 * flux_vector_[0];
                        modec_out_.width(18);

                        if (n_vector_[size_tot][2 * i] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[size_tot][2 * i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24* flux_vector_[size_tot];
                        }
                        modec_out_ << '\n';
                    }
                    if (print_mode_ == 1) {
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            if (n_vector_[j][2 * i] < cutoff || n_vector_[j][2 * i] < 0.0) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[j][2 * i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24* flux_vector_[j]; // 单位为Ci
                            }
                        }
                        modec_out_ << '\n';
                    }
                }
                modec_out_ << '\n';
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "Total: ";
                for (int j = 0; j <= size_tot; ++j) {
                    modec_out_.width(18);

                    modec_out_ << prod_neu[j]; // 单位为Ci

                }
                modec_out_ << '\n';

                modec_out_.close();

            }

            if (if_print_absorption_rate_ == 2) {
                output_filename_ = input_filename_ + ".NeuAbsRate";

                output_file_ = work_direc_ + output_filename_;
                modec_out_.open(output_file_);
                modec_out_ << "*********************************** Neutron Absorption Rate, s-1 ***********************************";
                modec_out_ << '\n';
                modec_out_ << '\n';

                modec_out_.setf(ios::left);
                int size = evolution_mode_.size();
                int size_tot = 0;
                for (int i = 0; i < size; ++i) {
                    size_tot += substep_[i];
                }

                vector<double> Cum_time;
                Cum_time.resize(size_tot + 1);

                int count = 0;
                for (int i = 0; i < size; ++i) {
                    for (int j = 0; j < substep_[i]; ++j) {
                        Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
                        count++;
                    }
                }


                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(18);
                modec_out_ << "Initial";

                if (print_mode_ == 0) {
                    modec_out_.width(10);
                    modec_out_.setf(ios::scientific | ios::uppercase);
                    modec_out_.precision(4);
                    modec_out_ << Cum_time[size_tot];
                    modec_out_.width(8);
                    modec_out_ << "s";
                } else {
                    for (int j = 1; j <= size_tot; ++j) {
                        modec_out_.width(10);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(4);
                        modec_out_ << Cum_time[j];
                        modec_out_.width(8);
                        modec_out_ << "s";
                    }
                }

                modec_out_ << '\n';
                modec_out_ << '\n';

                for (unsigned int i = 0; i < Nucl_size; ++i) {
                    modec_out_.width(8);
                    int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
                    modec_out_ << nucl_id;
                    int NZ = nucl_id / 10000;
                    int NA = (nucl_id - NZ * 10000) / 10;
                    stringstream ss;
                    ss << NA;
                    int NG = nucl_id - nucl_id / 10 * 10;
                    string Nucl_Name;
                    Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
                    if (NG == 1) {
                        Nucl_Name += "m";
                    }
                    modec_out_.width(8);
                    modec_out_ << Nucl_Name;
                    if (print_mode_ == 0) {
                        modec_out_.width(18);
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        modec_out_ << n_vector_[0][2 * i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24* flux_vector_[0];
                        modec_out_.width(18);
                        if (n_vector_[size_tot][2 * i] < cutoff) {
                            modec_out_ << 0.0;
                        } else {
                            modec_out_ << n_vector_[size_tot][2 * i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24* flux_vector_[size_tot];
                        }
                        modec_out_ << '\n';
                    }
                    if (print_mode_ == 1) {
                        modec_out_.setf(ios::scientific | ios::uppercase);
                        modec_out_.precision(9);
                        for (int j = 0; j <= size_tot; ++j) {
                            modec_out_.width(18);
                            if (n_vector_[j][2 * i] < cutoff) {
                                modec_out_ << 0.0;
                            } else {
                                modec_out_ << n_vector_[j][2 * i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24* flux_vector_[j];
                            }
                        }
                        modec_out_ << '\n';
                    }
                }

                modec_out_ << '\n';
                modec_out_.width(8);
                modec_out_ << "";
                modec_out_.width(8);
                modec_out_ << "Total: ";
                for (int j = 0; j <= size_tot; ++j) {
                    modec_out_.width(18);

                    modec_out_ << absorption_neu[j]; // 单位为Ci

                }
                modec_out_ << '\n';

                modec_out_.close();

            }

        }

    }
}
