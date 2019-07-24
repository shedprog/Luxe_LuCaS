

#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>

#include "LuxeTestGenerator.hh"


LuxeTestGenerator::LuxeTestGenerator() 
  : datf(0)  {

}

LuxeTestGenerator::~LuxeTestGenerator() {

  if (datf) delete datf;
}



void LuxeTestGenerator::AddEventFile (const std::string fname)
{
  if (!datf) {
    datf = new FileData();  
  }
  datf->AddFile(fname);
}


void LuxeTestGenerator::SetFileType (const std::string ftype, const int n_col, const int n_skip)
{
  if (!datf) {
    datf = new FileData();  
  }
  datf->SetFileType(ftype);
  datf->SetNColumns(n_col);
  datf->SetSkipLines(n_skip);
}


int LuxeTestGenerator::GetEventFromFile (std::vector < std::vector <double> > &ptcls)
{
  if (!datf) return -5;
  int nrdat = datf->GetData(ptcls);
  return nrdat;
}



void LuxeTestGenerator::SetFileList(const std::vector<std::string> fnlist)
{
  for (std::vector<std::string>::const_iterator itr = fnlist.begin(); itr != fnlist.end(); ++itr) {
    AddEventFile(*itr);  
  }
}



void LuxeTestGenerator::SetFileList(const std::string fnamelist)
{
  std::vector<std::string>  flist;
  ProcessList(fnamelist, flist);
  SetFileList(flist);
}



int LuxeTestGenerator::ProcessList(const std::string &fnamelist, std::vector<std::string> &flist)
{
  std::fstream  fdata;
  fdata.open(fnamelist, std::ios::in);
  if (!fdata.is_open()) {
    throw std::runtime_error(std::string("Error reding data from the file ") + fnamelist);
  }
  
  unsigned long lid = 0;
  while (!fdata.eof()) {
    std::string  ffname;
    fdata >> ffname;
    if (!fdata.fail()) { 
//       std::cout << "File name " << ffname << " is read from the list file" << std::endl;
      flist.push_back(ffname);
    }
    else if (fdata.eof()) { break; }
    else {
      std::cout << "ProcessList(..)  :  Error reading data from the file " << fnamelist 
                << ",  line: " << lid << ". Exit." << std::endl;
      fdata.close();          
      return -2;
    }
    ++lid;
  }
  
  fdata.close();

  return 0;
}  
  


//////////////////////////////////////////////////////////////////////////

FileData::FileData() : fcname(""), fid(-1), n_columns(17), skip_lines(0), freadfn(&FileData::ReadRecord), debugl(0)
{
}


FileData::FileData(const std::string fname) : fcname(""), fid(-1), n_columns(17), skip_lines(0), 
   freadfn(&FileData::ReadRecord), debugl(0)
{
  fname_list.push_back(fname);
}


FileData::FileData(const std::vector<std::string> &fnamelist) :
  fname_list(fnamelist), fcname(""), fid(-1), n_columns(17), skip_lines(0), 
  freadfn(&FileData::ReadRecord), debugl(0)
{
}


FileData::~FileData()
{
  if (fdata.is_open()) fdata.close();
}


int FileData::SetFileType(std::string ftype)
{
  if(!ftype.compare("hepevt")) {
    freadfn = &FileData::ReadRecord;
    n_columns = 18;
    skip_lines = 0;
  } else if (!ftype.compare("out")) { 
      freadfn = &FileData::ReadRecordOut;
      n_columns = 9;
      skip_lines = 9;
  } else {
    std::cout << "FileData::SetFileType: Error: File type " << ftype << " is not supported!\n";    
    return -1;
  }
  return 0;
}


void FileData::AddFile(const std::string fname)
{
  if (!fdata.is_open()) {
    fname_list.push_back(fname);
  } else {
    std::cout << "FileData::AddFile :  Cannot add file while reading data\n";
  }
}


int FileData::GetData(std::vector < std::vector <double> > &particles)
{
// It would be better probably to add exceptions for error handling not just return values
  int  ndata(0);  
  if (fdata.is_open() ) {
    ndata = (this->*freadfn)(particles);
    if (ndata >= 0 ) {
      return ndata;
    } else {  
      if ( fdata.eof() ) { fdata.close();  }
      else  { fdata.close(); return ndata; }
    }
  }

  if (!fdata.is_open()) {
    if ( ++fid < 0 || fname_list.size() > static_cast<unsigned int>(fid) ) {
      fdata.open(fname_list[fid], std::ios::in);
      fcname = fname_list[fid];
      if (debugl > 0) { std::cout << "FileData::GetData: opening file: " << fcname << std::endl; }
    } else {
      std::cout << "FileData::GetData(..) :  All data from " << fid << " file(s) have been read\n";
      fcname = "";
      --fid;
      return -1;
    }
  }
    
  if (!fdata.is_open()) {
    std::cerr << "FileData::GetData(..) :  Error open file " << fname_list[fid] << ". Exit.\n";
    fcname = "";
    return -2;
  }

  fdata.clear();
  std::string tmpstr;
  for (int ii = 0; ii < skip_lines; ++ii) {std::getline(fdata, tmpstr); }
  
  ndata = (this->*freadfn)(particles);
  if (ndata < 0) fdata.close();
      
  return ndata;
}


int FileData::ReadRecord(std::vector < std::vector <double> > &particles)
{
  std::vector <double> vdat(n_columns);
  int evid(0), npart(0);
   
  particles.clear(); 
  fdata >> evid >> npart;
  if (fdata.fail() || fdata.eof()) { 
    return -1; 
  }
  for (int ii = 0; ii < npart; ++ii) {
    for (int jj = 0; jj < n_columns; ++jj) {
      double xx;  
      fdata >> xx;
      if (!fdata.fail()) { vdat[jj] = xx; }
      else {
        std::cout << "FileData::ReadRecord(..)  :  Error reading data from the file " << fname_list[fid] 
                  << ",  event: " << evid << "  particle: " << ii << "  column: " << jj << std::endl;
        particles.clear();
        return -2;
      }
    }
    particles.push_back(vdat);
  }
  if (debugl > 1) {
    std::cout << "FileData::ReadRecord : Info: file: " << fcname << "  event: " << evid 
              << "  number of particles: " << npart << std::endl;
    for (auto itr = particles.begin(); itr != particles.end(); ++itr) {
      std::for_each(itr->begin(), itr->end(), [](const double x) {std::cout << x << "  ";} );
      std::cout << std::endl;
    }
  }
  return npart;
}


int FileData::ReadRecordOut(std::vector < std::vector <double> > &particles)
{
  std::vector <double> vdat(n_columns);
  int evid(0), npart(0);
   
  particles.clear(); 
  for (int jj = 0; jj < n_columns; ++jj) {
    double xx;  
    fdata >> xx;
//     std::cout << xx << "  ";
    if (!fdata.fail()) { vdat[jj] = xx; }
    else if (fdata.eof()) { return -1; }
    else {
      std::cout << "FileData::ReadRecordOut(..)  :  Error reading data from the file " << fname_list[fid] 
                << "  column: " << jj << std::endl;
      particles.clear();
      return -2;
    }
  }
  particles.push_back(vdat);
// std::cout << std::endl;

  if (debugl > 1) {
    std::cout << "FileData::ReadRecordOut : Info: file: " << fcname << "  event: " << evid << std::endl;
    for (auto itr = particles.begin(); itr != particles.end(); ++itr) {
      std::for_each(itr->begin(), itr->end(), [](const double x) {std::cout << x << "  ";} );
      std::cout << std::endl;
    }
  }
  return ++npart;
}


void FileData::ListFiles(void) const
{
  std::for_each(fname_list.begin(), fname_list.end(), [](const std::string ss) {std::cout << ss << std::endl;} );
}







