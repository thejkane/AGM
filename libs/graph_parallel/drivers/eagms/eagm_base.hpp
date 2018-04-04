// Copyright (C) 2018 Thejaka Amila Kanewala, Marcin Zalewski, Andrew Lumsdaine.

// Boost Software License - Version 1.0 - August 17th, 2003

// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:

// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

//  Authors: Thejaka Kanewala
//           Andrew Lumsdaine

#ifndef __EAGM_BASE__
#define __EAGM_BASE__
#include "../common/synthetic_generator.hpp"
#include "../common/graph_loader.hpp"
#include "../common/parser.hpp"
#include "../common/executor.hpp"

#include <boost/graph/util/work_stats.hpp>
#include <boost/graph/agm/model/eagm_buckets.hpp>
#include <boost/graph/agm/util/stat.hpp>
#include <boost/graph/agm/algorithms/sssp.hpp>
#include <boost/graph/agm/model/general_orderings.hpp>
#include <boost/graph/agm/runtime/ampp_runtime.hpp>

class agm_instance_params_base {
public:
  agm_instance_params_base(pf_execution_mode& _pf):pf_mode(_pf){}
  
  pf_execution_mode pf_mode;

  std::string get_pf_execution_mode() {
    if (pf_mode == agm_pf_preorder)
      return "pre-order";
    else if (pf_mode == agm_pf_postorder)
      return "post-order";
    else if (pf_mode == agm_pf_splitted)
      return "splitted";
    else {
      error("Invalid processing function mode. Available : preorder, postorder, splitted");
      assert(false);
    }            
  }

  void print() {
    std::cout << "processing function execution mode: " << get_pf_execution_mode()
              << std::endl;
  }
  
};

class agm_params_base {
protected:
  pf_execution_mode pf_mode = agm_pf_splitted;

  bool parse(int argc, char* argv[]) {
    for (int i=0; i < argc; ++i) {    
      std::string arg = argv[i];    
      if (arg == "--pf-exec-mode") {
        std::string execmode = argv[i+1];
        if (execmode == "preorder")
          pf_mode = agm_pf_preorder;
        else if (execmode == "postorder")
          pf_mode = agm_pf_postorder;
        else if (execmode == "splitted")
          pf_mode = agm_pf_splitted;
        else {
          error("Invalid processing function execution mode. Available (preorder, postorder, splitted)");
          assert(false);
        }
      }
    }
  }
};

#endif