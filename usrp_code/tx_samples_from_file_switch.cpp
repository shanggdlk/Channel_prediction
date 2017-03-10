#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/exception.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/thread.hpp>
#include <boost/chrono.hpp>
#include <iostream>
#include <fstream>
#include <complex>
#include <csignal>

namespace po = boost::program_options;

static bool stop_signal_called = false;
static bool freq_switch = false;
void sig_int_handler(int){
    freq_switch = true;
}

template<typename samp_type> void send_from_file(
    uhd::usrp::multi_usrp::sptr usrp,
    const std::string &cpu_format,
    const std::string &file,
    size_t samps_per_buff,
    double freq,
    double end_freq,
    double freq_step
){
    unsigned int num_channels = usrp->get_rx_num_channels();
    uhd::stream_args_t stream_args(cpu_format);
    for (size_t chan = 0; chan < num_channels; ++chan) {
        stream_args.channels.push_back(chan);
    }
    uhd::tx_streamer::sptr tx_stream = usrp->get_tx_stream(stream_args);
    uhd::time_spec_t timespec = usrp->get_time_now();


    uhd::tx_metadata_t md;
    md.start_of_burst = false;
    md.end_of_burst = true;
    md.has_time_spec = false;
    std::vector<samp_type> buff(samps_per_buff);
    std::vector<samp_type> zero_buff(samps_per_buff);
    std::vector<samp_type*> buffs(num_channels, &zero_buff.front());
    std::ifstream infile(file.c_str(), std::ifstream::binary);

    int is_end_of_file = 0;
    do {
        //loop until the entire file has been read
        while(!is_end_of_file) {

            infile.read((char*)&buff.front(), buff.size()*sizeof(samp_type));
            size_t num_tx_samps = infile.gcount()/sizeof(samp_type);

            is_end_of_file = infile.eof();
	        md.time_spec = timespec;

            buffs[0] = &buff.front();
            tx_stream->send(
                buffs, num_tx_samps, md
            );
        }
        infile.clear();
        infile.seekg(0, std::ios::beg);
        is_end_of_file = 0;
        if (freq_switch) {
            freq_switch = false;
            if (freq == end_freq) {
                stop_signal_called = true;
            }
            else {
                freq += freq_step;
                md.start_of_burst = false;
                std::cout << boost::format("Setting TX Freq: %f MHz...") % (freq/1e6) << std::endl;
                usrp->set_tx_freq(freq);
                std::cout << boost::format("Actual TX Freq: %f MHz...") % (usrp->get_tx_freq()/1e6) << std::endl << std::endl;
            }
        }
    } while (not stop_signal_called);

    infile.close();
}

int UHD_SAFE_MAIN(int argc, char *argv[]) {
    uhd::set_thread_priority_safe();

    //transmit variables to be set by po
    std::string tx_args, tx_file, tx_ant, tx_subdev, ref, type;
    double tx_rate, tx_freq, end_freq, freq_step, tx_gain, tx_bw;
    size_t spb;

    //setup the program options
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "help message")
        ("args", po::value<std::string>(&tx_args)->default_value(""), "uhd transmit device address args")
        ("file", po::value<std::string>(&tx_file)->default_value("usrp_samples.dat"), "name of the file to write binary samples to")
        ("type", po::value<std::string>(&type)->default_value("short"), "sample type in file: double, float, or short")
        ("spb", po::value<size_t>(&spb)->default_value(64000), "samples per buffer")
        ("rate", po::value<double>(&tx_rate), "rate of transmit outgoing samples")
        ("freq", po::value<double>(&tx_freq), "transmit RF center frequency in Hz")
        ("end_freq", po::value<double>(&end_freq), "RF end center frequency in Hz")
        ("freq_step", po::value<double>(&freq_step)->default_value(50000000), "RF frequency step in Hz")
        ("gain", po::value<double>(&tx_gain), "gain for the transmit RF chain")
        ("ant", po::value<std::string>(&tx_ant), "transmit antenna selection")
        ("subdev", po::value<std::string>(&tx_subdev), "transmit subdevice specification")
        ("bw", po::value<double>(&tx_bw), "analog transmit filter bandwidth in Hz")
        ("ref", po::value<std::string>(&ref)->default_value("internal"), "clock reference (internal, external, mimo)")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    //print the help message
    if (vm.count("help")){
        std::cout << boost::format("UHD RX samples to file %s") % desc << std::endl;
        return ~0;
    }

    //create a usrp device
    std::cout << std::endl;
    std::cout << boost::format("Creating the transmit usrp device with: %s...") % tx_args << std::endl;
    uhd::usrp::multi_usrp::sptr tx_usrp = uhd::usrp::multi_usrp::make(tx_args);

    //Lock mboard clocks
    if (ref == "mimo") {
        uhd::clock_config_t clock_config;
        clock_config.ref_source = uhd::clock_config_t::REF_MIMO;
        clock_config.pps_source = uhd::clock_config_t::PPS_MIMO;
        tx_usrp->set_clock_config(clock_config);
    }
    else if (ref == "external") {
        tx_usrp->set_clock_config(uhd::clock_config_t::external());
    }
    else if (ref == "internal") {
        tx_usrp->set_clock_config(uhd::clock_config_t::internal());
    }

    //always select the subdevice first, the channel mapping affects the other settings
    if (vm.count("subdev")) tx_usrp->set_tx_subdev_spec(tx_subdev);

    std::cout << boost::format("Using TX Device: %s") % tx_usrp->get_pp_string() << std::endl;

    const unsigned int tx_num_channels = tx_usrp->get_tx_num_channels();

    //set the transmit sample rate
    if (not vm.count("rate")){
        std::cerr << "Please specify the transmit sample rate with --tx-rate" << std::endl;
        return ~0;
    }
    std::cout << boost::format("Setting TX Rate: %f Msps...") % (tx_rate/1e6) << std::endl;
    for (unsigned int i = 0; i < tx_num_channels; ++i) {
        tx_usrp->set_tx_rate(tx_rate, i);
    }
    std::cout << boost::format("Actual TX Rate: %f Msps...") % (tx_usrp->get_tx_rate()/1e6) << std::endl << std::endl;

    //set the transmit center frequency
    if (not vm.count("freq")){
        std::cerr << "Please specify the transmit center frequency with --tx-freq" << std::endl;
        return ~0;
    }

    //uhd::time_spec_t tx_cmd_time = tx_usrp->get_time_now() + uhd::time_spec_t(0.1);
    //tx_usrp->set_command_time(tx_cmd_time);
    for (unsigned int i = 0; i < tx_num_channels; ++i) {
        std::cout << boost::format("Setting TX Freq: %f MHz...") % (tx_freq/1e6) << std::endl;
        uhd::tune_result_t r = tx_usrp->set_tx_freq(tx_freq, i);
        std::cout << r.to_pp_string() << std::endl;
        std::cout << boost::format("Actual TX Freq: %f MHz...") % (tx_usrp->get_tx_freq()/1e6) << std::endl << std::endl;
    }
    //tx_usrp->clear_command_time();

    //set the transmit rf gain
    if (vm.count("gain")){
        for (unsigned int i = 0; i < tx_num_channels; ++i) {
            std::cout << boost::format("Setting TX Gain: %f dB...") % tx_gain << std::endl;
            tx_usrp->set_tx_gain(tx_gain, i);
            std::cout << boost::format("Actual TX Gain: %f dB...") % tx_usrp->get_tx_gain() << std::endl << std::endl;
        }
    }

    //set the analog frontend filter bandwidth
    if (vm.count("bw")){
        for (unsigned int i = 0; i < tx_num_channels; ++i) {
            std::cout << boost::format("Setting TX Bandwidth: %f MHz...") % tx_bw << std::endl;
            tx_usrp->set_tx_bandwidth(tx_bw, i);
            std::cout << boost::format("Actual TX Bandwidth: %f MHz...") % tx_usrp->get_tx_bandwidth() << std::endl << std::endl;
        }
    }

    //set the antenna
    if (vm.count("ant")) {
        for (unsigned int i = 0; i < tx_num_channels; ++i) {
            tx_usrp->set_tx_antenna(tx_ant, i);
        }
    }

    boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    //Check Ref and LO Lock detect
    std::vector<std::string> tx_sensor_names;
    tx_sensor_names = tx_usrp->get_tx_sensor_names(0);
    if (std::find(tx_sensor_names.begin(), tx_sensor_names.end(), "lo_locked") != tx_sensor_names.end()) {
        uhd::sensor_value_t lo_locked = tx_usrp->get_tx_sensor("lo_locked",0);
        std::cout << boost::format("Checking TX: %s ...") % lo_locked.to_pp_string() << std::endl;
        UHD_ASSERT_THROW(lo_locked.to_bool());
    }

    tx_sensor_names = tx_usrp->get_mboard_sensor_names(0);
    if ((ref == "mimo") and (std::find(tx_sensor_names.begin(), tx_sensor_names.end(), "mimo_locked") != tx_sensor_names.end())) {
        uhd::sensor_value_t mimo_locked = tx_usrp->get_mboard_sensor("mimo_locked",0);
        std::cout << boost::format("Checking TX: %s ...") % mimo_locked.to_pp_string() << std::endl;
        UHD_ASSERT_THROW(mimo_locked.to_bool());
    }
    if ((ref == "external") and (std::find(tx_sensor_names.begin(), tx_sensor_names.end(), "ref_locked") != tx_sensor_names.end())) {
        uhd::sensor_value_t ref_locked = tx_usrp->get_mboard_sensor("ref_locked",0);
        std::cout << boost::format("Checking TX: %s ...") % ref_locked.to_pp_string() << std::endl;
        UHD_ASSERT_THROW(ref_locked.to_bool());
    }

    std::signal(SIGINT, &sig_int_handler);
    std::cout << "Press Ctrl + C to stop streaming..." << std::endl;

    //send from file
    if (type == "double") send_from_file<std::complex<double> >(tx_usrp, "fc64", tx_file, spb, tx_freq, end_freq, freq_step);
    else if (type == "float") send_from_file<std::complex<float> >(tx_usrp, "fc32", tx_file, spb, tx_freq, end_freq, freq_step);
    else if (type == "short") send_from_file<std::complex<short> >(tx_usrp, "sc16", tx_file, spb, tx_freq, end_freq, freq_step);
    else {
      throw std::runtime_error("Unknown type " + type);
    }

    //finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;
    return 0;
}
