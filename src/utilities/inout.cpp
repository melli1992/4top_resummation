#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include "parameters.h"
#include <iostream>
#include <fstream>
#include <iterator>
using namespace std;

//////////////////////////////////////////////////////////////////
///
/// modification of input parameters
/// default file to read input from: input.cfg
///
//////////////////////////////////////////////////////////////////
// A helper function to simplify the main part.
template<class T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
    return os;
}


int configure(int ac, char* av[], string configfile, bool runupdate = true)
{
    try {

        double comTeV; // sqrt(S) in TeV
        string config_file = configfile;
		// Declare a group of options that will be
        // allowed only on command line
        po::options_description generic("Generic options");
        generic.add_options()
            ("help,h", "produce help message")
            ("config,c", po::value<string>(&config_file)->default_value(configfile), "name of a file of a configuration.")
            ;

        // Declare a group of options that will be
        // allowed both on command line and in
        // config file
        po::options_description config("Configuration");
        config.add_options()
            ("include-path,I",
                 po::value< vector<string> >()->composing(),
                 "include path")
            ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        po::options_description hidden("Hidden options");
        hidden.add_options()
            //PDFstuf
            ("setname", po::value<string>(&setname)->default_value("PDF4LHC15_nnlo_100"), "Name of PDFset")
            ("usemember", po::value<int>(&use_member)->default_value(0),"member of PDFset to use")
            ("fitPDF", po::value<bool>(&fitPDF)->default_value(false),"use own fits to PDF grid")
            ("method", po::value<string>(&method)->default_value("res"), "method of integration")
            //process
            ("sqrtS", po::value<double>(&comTeV)->default_value(13.), "center-of-mass energy [TeV]")
            //("observable", po::value<string>(&observable)->default_value("Qinv"),"observable")
            //resummation
            ("inceuler", po::value<double>(&INCEULER)->default_value(1), "resum euler constant in pQCD")
            ("expansion", po::value<bool>(&expansion)->default_value(false), "expanded result up to NLO")
            ("includeS1", po::value<bool>(&include_S1)->default_value(false), "include S1")
            ("includeC1", po::value<bool>(&include_C1)->default_value(false), "include C1")
            ("qqchan", po::value<bool>(&include_qqbar)->default_value(true), "include qqbar")
            ("ggchan", po::value<bool>(&include_gg)->default_value(true), "include gg")
            //scales
            ("muR", po::value<double>(&muR)->default_value(500.), "renormalization scale [GeV]")
            ("muF", po::value<double>(&muF)->default_value(500.), "factorization scale [GeV]")
            ("Q", po::value<double>(&Q)->default_value(500.), "hard scale [GeV]")
            //masses
            ("mt", po::value<double>(&mt)->default_value(172.5), "Top quark mass [GeV]")
            //integration
            ("phiMP", po::value<double>(&phiMP)->default_value(3./4.*M_PI), "phiMP for inverse Mellin transform")
            ("CMP", po::value<double>(&CMP)->default_value(2.1), "CMP for inverse Mellin transform")
            ;


        po::options_description cmdline_options;
        cmdline_options.add(generic).add(config).add(hidden);

        po::options_description config_file_options;
        config_file_options.add(config).add(hidden);

        po::options_description visible("Allowed options");
        visible.add(generic).add(config);

                po::positional_options_description p;
        p.add("input-file", -1);

        po::variables_map vm;
        store(po::command_line_parser(ac, av).
              options(cmdline_options).positional(p).run(), vm);
        notify(vm);

        ifstream ifs(config_file.c_str());
        if (!ifs)
        {
            cout << "Cannot open config file: " << config_file << "\n";
            exit(0);
        }
        else
        {
            store(parse_config_file(ifs, config_file_options), vm);
            notify(vm);
        }

        if (vm.count("help")) {
			  cout << endl
				   << "fill in the input.cfg file" << endl;
			  exit(0);
        }
        if (vm.count("sqrtS")){
            S = comTeV*1000.;
            S2 = pow(S,2);
        }
        if( vm.count("mt")){mt2 = pow(mt,2);}
        if (runupdate){update_defaults();}

    }
    catch(exception& e)
    {
        cout << e.what() << "\n";
        return 1;
    }
    return 0;
}
     
