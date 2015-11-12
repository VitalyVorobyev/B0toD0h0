#include "filter.h"

int main(int argc, char** argv){
  int type_in = 0;
  int type_out = 0;
  int stream = 0;
  bool line_flag = false;
  if(argc == 2){
    sscanf(argv[1],"%d",&type_in);
    switch(type_in){
    case 11: type_out = 1; break;
    case 12: type_out = 2; break;
    case 13: type_out = 4; break;
    case 15: type_out = 5; break;
    case 101: type_out = 10; break;
    }
    type_out == type_in;
  } else if(argc == 3){
    sscanf(argv[1],"%d",&type_in);
    sscanf(argv[2],"%d",&stream);
    type_out == type_in;
  } else if(argc == 4){
    sscanf(argv[1],"%d",&type_in);
    sscanf(argv[2],"%d",&stream);
    sscanf(argv[3],"%d",&type_out);
  } else if(argc == 5){
    sscanf(argv[1],"%d",&type_in);
    sscanf(argv[2],"%d",&stream);
    sscanf(argv[3],"%d",&type_out);
    line_flag = true;
  }
  if(line_flag){
    vector<string> angle_vec;
    angle_vec.push_back(string("0"));
    angle_vec.push_back(string("15"));
    angle_vec.push_back(string("22_5"));
    angle_vec.push_back(string("30"));
    angle_vec.push_back(string("45"));
    angle_vec.push_back(string("60"));
    angle_vec.push_back(string("67_5"));
    angle_vec.push_back(string("75"));
    angle_vec.push_back(string("90"));
    angle_vec.push_back(string("105"));
    angle_vec.push_back(string("112_5"));
    angle_vec.push_back(string("120"));
    angle_vec.push_back(string("135"));
    angle_vec.push_back(string("150"));
    angle_vec.push_back(string("157_5"));
    angle_vec.push_back(string("165"));
    for(int i=0; i<angle_vec.size(); i++){
      filter fil(type_in,type_out,stream,angle_vec[i]);
      fil.MakeNTuples();
    }
    return 0;
  } else{
    filter fil(type_in,type_out,stream);
    fil.MakeNTuples();
    return 0;
  }

  return 0;
}
