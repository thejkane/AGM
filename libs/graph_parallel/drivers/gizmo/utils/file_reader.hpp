#ifndef __AGM_GIZMO_FILE_READER__
#define __AGM_GIZMO_FILE_READER__

#define KEY_BANDWIDTH "bandwidth"
#define KEY_RANKS "ranks"
#define KEY_THREADS "threads"
#define KEY_LATENCY_PROFILE "latency-profile"
#define KEY_SEND_OVERHEAD "send-overhead"
#define KEY_RECV_OVERHEAD "receive-overhead"
#define KEY_LATENCY_PROFILE_ENABLED "latency-profile.enabled"
#define KEY_LATENCY_VALUE "latency-value"
#define KEY_SELF_SEND_ENABLED "self-send-enabled"

class GizmoConfigException : public std::exception {
private: 
  std::string err;
public:
  GizmoConfigException(const char* e) : err(e){
    std::cout << e << std::endl;
  }

  GizmoConfigException(std::string e) : err(e){
    std::cout << e << std::endl;
  }

  const char * what () const throw () {
    return err.c_str();
  }  
};


enum machine_model {
  distributed_real,
  pram,
  ram,
  invalid
};


class gizmo_config {
private:
  std::map<std::string, std::string> configs;
  machine_model machine;
public:
  void insert(std::string k, std::string v) {
    configs.insert(std::make_pair(k, v));
  }

  void set_machine_model(machine_model model) {
    machine = model;
  }

  std::string get(std::string k) const {
    auto iter = configs.find(k);
    if (iter != configs.end()) {
      return iter->second;
    } else {
      std::string serr = "Key : " + k + " not found.";
      throw GizmoConfigException(serr);
    }
  }

  bool get_bool(std::string k) const {
    return boost::lexical_cast<bool>(get(k));
  }

  double get_double(std::string k) const {
    return std::stod(get(k));
  }

  int get_int(std::string k) const {
    int i = std::stoi(get(k));
    return i;
  }

  machine_model get_machine_model() {
    return machine;
  }

  void print() const {
    std::cout << "================ Gizmo Configurations =====================" << std::endl;
    std::map<std::string, std::string>::const_iterator
      itebegin = configs.begin();
    for (; itebegin != configs.end(); ++itebegin) {
      std::cout << (*itebegin).first << " : " 
		<< (*itebegin).second << std::endl;
    }

    std::cout << "===========================================================" << std::endl;
  }
};


class config_file_reader {
private:
  std::string file;
  machine_model machine;
  
  machine_model get_machine_model(std::string _m) {
    if (_m == "distributed_real")
      return distributed_real;
    else if (_m == "pram")
      return pram;
    else if (_m == "ram")
      return ram;
    else {
      std::cout << "[ERROR] Invalid machine type. Available types are :"
		<< " distributed_real, pram and ram" << std::endl;
      return invalid;
    }
  }

public:
  config_file_reader(const char* f) : file(f) {}
  config_file_reader() : file("") {}

  bool parse(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];

      if (arg == "--gizmo-file") {
	file = boost::lexical_cast<std::string>( argv[i+1]);
      }

      if (arg == "--machine") {
	std::string smmodel = argv[i+1];
	machine = get_machine_model(smmodel);
	if (machine == invalid)
	  return false;
      }
    }

    if (file == "") {
      std::cout << "[ERROR] Gizmo configuration file not found. Use --gizmo-file to specify the configuration file." << std::endl;
      exit(-1);
    }

    return true;
  }

  void read_configurations(gizmo_config& configs) {

    configs.set_machine_model(machine);
    std::ifstream ffile(file);

    if (!ffile) {
      std::string serr = "[ERROR] Error reading the file : " + file; 
      ffile.close();
      throw GizmoConfigException(serr);
    }
    
    std::string line;
    
    while( std::getline(ffile, line)) {
      std::istringstream is_line(line);
      std::string key;
      if( std::getline(is_line, key, '=')) {
	if (key[0]=='#')
	  continue;

	std::string value;
	if( std::getline(is_line, value) ) {
	  //	  std::cout << "inserting : " << key << ", " << value << std::endl;
	  configs.insert(key, value);
	} else {
	  std::string serr = "[ERROR] Configuration for key : " + key + ", is invalid."; 
	  ffile.close();
	  throw GizmoConfigException(serr);
	}
      }
    }
    
    ffile.close();
  }
};

#endif
