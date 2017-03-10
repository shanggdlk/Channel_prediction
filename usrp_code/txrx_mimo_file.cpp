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
void sig_int_handler(int){
    stop_signal_called = true;
}

template<typename samp_type> void send_from_file(
    uhd::usrp::multi_usrp::sptr usrp,
    const std::string &cpu_format,
    const std::string &file,
    size_t samps_per_buff,
    double tx_offset,
    double tx_duration
){
    unsigned int num_channels = usrp->get_rx_num_channels();
    uhd::stream_args_t stream_args(cpu_format);
    for (size_t chan = 0; chan < num_channels; ++chan) {
        stream_args.channels.push_back(chan);
    }
    uhd::tx_streamer::sptr tx_stream = usrp->get_tx_stream(stream_args);
    uhd::time_spec_t timespec = usrp->get_time_now();
    uhd::time_spec_t start_offset(tx_offset);
    uhd::time_spec_t duration(tx_duration);
    timespec += start_offset;


    uhd::tx_metadata_t md;
    md.start_of_burst = true;
    md.end_of_burst = true;
    md.has_time_spec = (tx_duration != 0);
    std::vector<samp_type> buff(samps_per_buff);
    std::vector<samp_type> zero_buff(samps_per_buff);
    std::vector<samp_type*> buffs(num_channels, &zero_buff.front());
    std::ifstream infile(file.c_str(), std::ifstream::binary);

    int is_end_of_file = 0;
    while (not stop_signal_called) {
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

            for (unsigned int i = 1; i < num_channels; ++i) {
                md.time_spec += 0.001;
                buffs[i - 1] = &zero_buff.front();
                buffs[i] = &buff.front();
                tx_stream->send(
                    buffs, num_tx_samps, md
                );
            }
        }
        infile.clear();
        infile.seekg(0, std::ios::beg);
        is_end_of_file = 0;
        timespec += duration;
    }

    infile.close();
}

template<typename samp_type> void recv_to_file(
    uhd::usrp::multi_usrp::sptr usrp,
    const std::string &cpu_format,
    const std::string &file,
    size_t samps_per_buff,
    double seconds_in_future,
    int num_requested_samples
){
    int num_total_samps = 0;
    unsigned int num_channels = usrp->get_rx_num_channels();
    uhd::stream_args_t stream_args(cpu_format, "sc16");
    for (size_t chan = 0; chan < num_channels; ++chan) {
        stream_args.channels.push_back(chan);
    }
    uhd::rx_streamer::sptr rx_stream = usrp->get_rx_stream(stream_args);
    uhd::rx_metadata_t md;
    std::vector<samp_type *>buff(num_channels);
    std::ofstream outfiles[num_channels];

    for (unsigned int i = 0; i < num_channels; ++i) {
      buff[i] = new samp_type[samps_per_buff];
      std::string rx_file = file + "_" + boost::lexical_cast<std::string>(i) + ".dat";
      outfiles[i].open(rx_file.c_str(), std::ofstream::binary);
      std::cout << boost::format("Channel %u: Writing to file %s...\n") % i % rx_file;
    }

    bool overflow_message = true;

    uhd::stream_cmd_t stream_cmd((num_requested_samples == 0)?
        uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS:
        uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE
    );
    //uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);
    stream_cmd.num_samps = num_requested_samples;
    //stream_cmd.num_samps = 0;
    stream_cmd.stream_now = false;
    stream_cmd.time_spec = usrp->get_time_now() + uhd::time_spec_t(seconds_in_future);
    usrp->issue_stream_cmd(stream_cmd);

    while(not stop_signal_called and (num_requested_samples > num_total_samps or num_requested_samples == 0)){
            size_t num_rx_samps = rx_stream->recv(
            buff, samps_per_buff, md,
            seconds_in_future + 4.0
        );
        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) {
            std::cout << boost::format("Timeout while streaming") << std::endl;
            break;
        }
        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW){
            if (overflow_message){
                overflow_message = false;
                std::cerr << boost::format(
                    "Got an overflow indication. Please consider the following:\n"
                    "  Your write medium must sustain a rate of %fMB/s.\n"
                    "  Dropped samples will not be written to the file.\n"
                    "  Please modify this example for your purposes.\n"
                    "  This message will not appear again.\n"
                ) % (usrp->get_rx_rate()*sizeof(samp_type)/1e6);
            }
            continue;
        }
        if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE){
            throw std::runtime_error(str(boost::format(
                "Unexpected error code 0x%x"
            ) % md.error_code));
        }
        num_total_samps += num_rx_samps;
        //outfile.write((const char*)&buff.front(), num_rx_samps*sizeof(samp_type));
	    for (unsigned int i = 0; i < num_channels; ++i) {
	        outfiles[i].write((const char*)(buff[i]), num_rx_samps*sizeof(samp_type));
	    }
    }
    for (unsigned int i = 0; i < num_channels; ++i) {
      outfiles[i].close();
    }
    //outfile.close();
    stop_signal_called = true;
}

void receive_worker(
    uhd::usrp::multi_usrp::sptr usrp,
    const std::string &type,
    const std::string &file,
    size_t samps_per_buff,
    double seconds_in_future,
    int num_requested_samples
){
    if (type == "double") recv_to_file<std::complex<double> >(usrp, "fc64", file, samps_per_buff, seconds_in_future, num_requested_samples);
    else if (type == "float") recv_to_file<std::complex<float> >(usrp, "fc32", file, samps_per_buff, seconds_in_future, num_requested_samples);
    else if (type == "short") recv_to_file<std::complex<short> >(usrp, "sc16", file, samps_per_buff, seconds_in_future, num_requested_samples);
    else throw std::runtime_error("Unknown type " + type);
}


int UHD_SAFE_MAIN(int argc, char *argv[]) {
    uhd::set_thread_priority_safe();

    //transmit variables to be set by po
    std::string tx_args, tx_file, tx_ant, tx_subdev, ref;
    double tx_rate, tx_freq, tx_gain, tx_bw, tx_duration, tx_offset;

    //receive variables to be set by po
    std::string rx_args, rx_file, type, rx_ant, rx_subdev;
    double rx_rate, rx_freq, rx_gain, rx_bw;
    size_t spb, total_num_samps;

    float settling;

    //setup the program options
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "help message")
        ("tx-args", po::value<std::string>(&tx_args)->default_value(""), "uhd transmit device address args")
        ("rx-args", po::value<std::string>(&rx_args)->default_value(""), "uhd receive device address args")
        ("tx-file", po::value<std::string>(&tx_file)->default_value("usrp_samples.dat"), "name of the file to write binary samples to")
        ("rx-file", po::value<std::string>(&rx_file)->default_value("usrp_samples2.dat"), "name of the file to write binary samples to")
        ("type", po::value<std::string>(&type)->default_value("short"), "sample type in file: double, float, or short")
        ("settling", po::value<float>(&settling)->default_value(float(0.2)), "settling time (seconds) before receiving")
        ("spb", po::value<size_t>(&spb)->default_value(64000), "samples per buffer")
        ("tx-rate", po::value<double>(&tx_rate), "rate of transmit outgoing samples")
        ("rx-rate", po::value<double>(&rx_rate), "rate of receive incoming samples")
        ("tx-freq", po::value<double>(&tx_freq), "transmit RF center frequency in Hz")
        ("rx-freq", po::value<double>(&rx_freq), "receive RF center frequency in Hz")
        ("tx-gain", po::value<double>(&tx_gain), "gain for the transmit RF chain")
        ("rx-gain", po::value<double>(&rx_gain), "gain for the receive RF chain")
        ("tx-ant", po::value<std::string>(&tx_ant), "transmit antenna selection")
        ("rx-ant", po::value<std::string>(&rx_ant), "receive antenna selection")
        ("tx-subdev", po::value<std::string>(&tx_subdev), "transmit subdevice specification")
        ("rx-subdev", po::value<std::string>(&rx_subdev), "receive subdevice specification")
        ("tx-bw", po::value<double>(&tx_bw), "analog transmit filter bandwidth in Hz")
        ("rx-bw", po::value<double>(&rx_bw), "analog receive filter bandwidth in Hz")
        ("ref", po::value<std::string>(&ref)->default_value("internal"), "clock reference (internal, external, mimo)")
        ("nsamps", po::value<size_t>(&total_num_samps)->default_value(0), "total number of samples to receive")
        ("tx-duration", po::value<double>(&tx_duration)->default_value(0), "the time separation of each packet to transmit")
        ("tx-offset", po::value<double>(&tx_offset)->default_value(1), "the start of transmitter after the receiver started")
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
    std::cout << std::endl;
    std::cout << boost::format("Creating the receive usrp device with: %s...") % rx_args << std::endl;
    uhd::usrp::multi_usrp::sptr rx_usrp = uhd::usrp::multi_usrp::make(rx_args);

    //Lock mboard clocks
    if (ref == "mimo") {
        uhd::clock_config_t clock_config;
        clock_config.ref_source = uhd::clock_config_t::REF_MIMO;
        clock_config.pps_source = uhd::clock_config_t::PPS_MIMO;
        tx_usrp->set_clock_config(clock_config);
        rx_usrp->set_clock_config(clock_config);
    }
    else if (ref == "external") {
        tx_usrp->set_clock_config(uhd::clock_config_t::external());
        rx_usrp->set_clock_config(uhd::clock_config_t::external());
       // tx_usrp->set_time_unknown_pps(uhd::time_spec_t(0.0));
        rx_usrp->set_time_unknown_pps(uhd::time_spec_t(0.0));
        boost::this_thread::sleep(boost::posix_time::seconds(1)); //allow for some setup time
    }
    else if (ref == "internal") {
        tx_usrp->set_clock_config(uhd::clock_config_t::internal());
        rx_usrp->set_clock_config(uhd::clock_config_t::internal());
    }

    //always select the subdevice first, the channel mapping affects the other settings
    if (vm.count("tx-subdev")) tx_usrp->set_tx_subdev_spec(tx_subdev);
    if (vm.count("rx-subdev")) rx_usrp->set_rx_subdev_spec(rx_subdev);

    std::cout << boost::format("Using TX Device: %s") % tx_usrp->get_pp_string() << std::endl;
    std::cout << boost::format("Using RX Device: %s") % rx_usrp->get_pp_string() << std::endl;

    const unsigned int tx_num_channels = tx_usrp->get_tx_num_channels();
    const unsigned int rx_num_channels = rx_usrp->get_rx_num_channels();

    //set the transmit sample rate
    if (not vm.count("tx-rate")){
        std::cerr << "Please specify the transmit sample rate with --tx-rate" << std::endl;
        return ~0;
    }
    std::cout << boost::format("Setting TX Rate: %f Msps...") % (tx_rate/1e6) << std::endl;
    for (unsigned int i = 0; i < tx_num_channels; ++i) {
        tx_usrp->set_tx_rate(tx_rate, i);
    }
    std::cout << boost::format("Actual TX Rate: %f Msps...") % (tx_usrp->get_tx_rate()/1e6) << std::endl << std::endl;

    //set the receive sample rate
    if (not vm.count("rx-rate")){
        std::cerr << "Please specify the sample rate with --rx-rate" << std::endl;
        return ~0;
    }
    std::cout << boost::format("Setting RX Rate: %f Msps...") % (rx_rate/1e6) << std::endl;
    for (unsigned int i = 0; i < rx_num_channels; ++i) {
        rx_usrp->set_rx_rate(rx_rate, i);
    }
    std::cout << boost::format("Actual RX Rate: %f Msps...") % (rx_usrp->get_rx_rate()/1e6) << std::endl << std::endl;

    //set the transmit center frequency
    if (not vm.count("tx-freq")){
        std::cerr << "Please specify the transmit center frequency with --tx-freq" << std::endl;
        return ~0;
    }

    uhd::time_spec_t tx_cmd_time = tx_usrp->get_time_now() + uhd::time_spec_t(0.1);
    tx_usrp->set_command_time(tx_cmd_time);
    for (unsigned int i = 0; i < tx_num_channels; ++i) {
        std::cout << boost::format("Setting TX Freq: %f MHz...") % (tx_freq/1e6) << std::endl;
        uhd::tune_result_t r = tx_usrp->set_tx_freq(tx_freq, i);
        std::cout << r.to_pp_string() << std::endl;
        std::cout << boost::format("Actual TX Freq: %f MHz...") % (tx_usrp->get_tx_freq()/1e6) << std::endl << std::endl;
    }
    tx_usrp->clear_command_time();

    //set the receive center frequency
    if (not vm.count("rx-freq")){
        std::cerr << "Please specify the receive center frequency with --rx-freq" << std::endl;
        return ~0;
    }

    uhd::time_spec_t rx_cmd_time = rx_usrp->get_time_now() + uhd::time_spec_t(0.1);
    rx_usrp->set_command_time(rx_cmd_time);
    for (unsigned int i = 0; i < rx_num_channels; ++i) {
        std::cout << boost::format("Setting RX Freq: %f MHz...") % (rx_freq/1e6) << std::endl;
        uhd::tune_result_t r = rx_usrp->set_rx_freq(rx_freq, i);
        std::cout << r.to_pp_string() << std::endl;
        std::cout << boost::format("Actual RX Freq: %f MHz...") % (rx_usrp->get_rx_freq()/1e6) << std::endl << std::endl;
    }
    rx_usrp->clear_command_time();

    //set the transmit rf gain
    if (vm.count("tx-gain")){
        for (unsigned int i = 0; i < tx_num_channels; ++i) {
            std::cout << boost::format("Setting TX Gain: %f dB...") % tx_gain << std::endl;
            tx_usrp->set_tx_gain(tx_gain, i);
            std::cout << boost::format("Actual TX Gain: %f dB...") % tx_usrp->get_tx_gain() << std::endl << std::endl;
        }
    }

    //set the receive rf gain
    if (vm.count("rx-gain")){
        for (unsigned int i = 0; i < rx_num_channels; ++i) {
            std::cout << boost::format("Setting RX Gain: %f dB...") % rx_gain << std::endl;
            rx_usrp->set_rx_gain(rx_gain, i);
            std::cout << boost::format("Actual RX Gain: %f dB...") % rx_usrp->get_rx_gain() << std::endl << std::endl;
        }
    }

    //set the analog frontend filter bandwidth
    if (vm.count("tx-bw")){
        for (unsigned int i = 0; i < tx_num_channels; ++i) {
            std::cout << boost::format("Setting TX Bandwidth: %f MHz...") % tx_bw << std::endl;
            tx_usrp->set_tx_bandwidth(tx_bw, i);
            std::cout << boost::format("Actual TX Bandwidth: %f MHz...") % tx_usrp->get_tx_bandwidth() << std::endl << std::endl;
        }
    }

    //set the receive analog frontend filter bandwidth
    if (vm.count("rx-bw")){
        for (unsigned int i = 0; i < rx_num_channels; ++i) {
            std::cout << boost::format("Setting RX Bandwidth: %f MHz...") % (rx_bw/1e6) << std::endl;
            rx_usrp->set_rx_bandwidth(rx_bw, i);
            std::cout << boost::format("Actual RX Bandwidth: %f MHz...") % (rx_usrp->get_rx_bandwidth()/1e6) << std::endl << std::endl;
        }
    }

    //set the antenna
    if (vm.count("tx-ant")) {
        for (unsigned int i = 0; i < tx_num_channels; ++i) {
            tx_usrp->set_tx_antenna(tx_ant, i);
        }
    }
    //set the receive antenna
    if (vm.count("rx-ant")) {
        for (unsigned int i = 0; i < rx_num_channels; ++i) {
            rx_usrp->set_rx_antenna(rx_ant, i);
        }
    }

    //Check Ref and LO Lock detect
    boost::this_thread::sleep_for(boost::chrono::milliseconds(100));
    std::vector<std::string> tx_sensor_names, rx_sensor_names;
    tx_sensor_names = tx_usrp->get_tx_sensor_names(0);
    if (std::find(tx_sensor_names.begin(), tx_sensor_names.end(), "lo_locked") != tx_sensor_names.end()) {
        uhd::sensor_value_t lo_locked = tx_usrp->get_tx_sensor("lo_locked",0);
        std::cout << boost::format("Checking TX: %s ...") % lo_locked.to_pp_string() << std::endl;
        UHD_ASSERT_THROW(lo_locked.to_bool());
    }
    rx_sensor_names = rx_usrp->get_rx_sensor_names(0);
    if (std::find(rx_sensor_names.begin(), rx_sensor_names.end(), "lo_locked") != rx_sensor_names.end()) {
        uhd::sensor_value_t lo_locked = rx_usrp->get_rx_sensor("lo_locked",0);
        std::cout << boost::format("Checking RX: %s ...") % lo_locked.to_pp_string() << std::endl;
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

    rx_sensor_names = rx_usrp->get_mboard_sensor_names(0);
    if ((ref == "mimo") and (std::find(rx_sensor_names.begin(), rx_sensor_names.end(), "mimo_locked") != rx_sensor_names.end())) {
        uhd::sensor_value_t mimo_locked = rx_usrp->get_mboard_sensor("mimo_locked",0);
        std::cout << boost::format("Checking RX: %s ...") % mimo_locked.to_pp_string() << std::endl;
        UHD_ASSERT_THROW(mimo_locked.to_bool());
    }
    if ((ref == "external") and (std::find(rx_sensor_names.begin(), rx_sensor_names.end(), "ref_locked") != rx_sensor_names.end())) {
        uhd::sensor_value_t ref_locked = rx_usrp->get_mboard_sensor("ref_locked",0);
        std::cout << boost::format("Checking RX: %s ...") % ref_locked.to_pp_string() << std::endl;
        UHD_ASSERT_THROW(ref_locked.to_bool());
    }

    std::signal(SIGINT, &sig_int_handler);
    std::cout << "Press Ctrl + C to stop streaming..." << std::endl;

    // start receive worker thread
    boost::thread_group receive_thread;
    try {
        receive_thread.create_thread(boost::bind(&receive_worker, rx_usrp, type, rx_file, spb, settling, total_num_samps));
    } catch(...) {
        stop_signal_called = true;
        receive_thread.join_all();
    }

    //send from file
    if (type == "double") send_from_file<std::complex<double> >(tx_usrp, "fc64", tx_file, spb, tx_offset, tx_duration);
    else if (type == "float") send_from_file<std::complex<float> >(tx_usrp, "fc32", tx_file, spb, tx_offset, tx_duration);
    else if (type == "short") send_from_file<std::complex<short> >(tx_usrp, "sc16", tx_file, spb, tx_offset, tx_duration);
    else {
      stop_signal_called = true;
      receive_thread.join_all();
      throw std::runtime_error("Unknown type " + type);
    }
    //clean up transmit worker
    stop_signal_called = true;
    receive_thread.join_all();
    //finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;
    return 0;
}
