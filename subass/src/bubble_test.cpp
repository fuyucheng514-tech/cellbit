#include "bubbles.hpp"
#include <iostream>
int main(int argc,char**argv){
  if(argc!=4){std::cerr<<"usage: bubble-test in.bam ref.fa out.fa\n";return 2;}
  try{auto r=subass::make_bubbles(argv[1],argv[2],argv[3],1);std::cout<<"contigs="<<r.mean_coverage.size()<<" error="<<r.mean_alignment_error<<"\n";return 0;}
  catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}
}
