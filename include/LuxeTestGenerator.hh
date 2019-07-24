
#ifndef __LUXETESTGENERATOR_H
#define __LUXETESTGENERATOR_H

#include <fstream>
#include <string>
#include <vector>


class FileData;  // Described heareafter

class LuxeTestGenerator {

public:  

  /// Constructor  
 // LuxeTestGenerator (Laser &laser_, BaseEBeam &ebeam, int rndseed_);
  LuxeTestGenerator ();
  
  /// Destructor
  ~LuxeTestGenerator ();

  /// do it

  
  void AddEventFile (const std::string fname);
  void SetFileType (const std::string ftype, const int n_col = 17, const int n_skip = 0);
  void SetFileList(const std::vector<std::string> fnlist);
  void SetFileList(const std::string fnamelist);
  int GetEventFromFile (std::vector < std::vector <double> > &ptcls);

protected:  
  int ProcessList(const std::string &fnamelist, std::vector<std::string> &flist);
  
protected:

  FileData  *datf;
};  


class FileData
{
public:    
  FileData();
  FileData(const std::string fname);
  FileData(const std::vector<std::string> &fnamelist);
  virtual ~FileData();
  void AddFile(const std::string fname);
  void SetNColumns (const int n) { n_columns = n; }
  void SetSkipLines (const int n) { skip_lines = n; }
  void SetDebug (const int n) { debugl = n; }

  int GetNColumns () const { return n_columns; }
  int GetData(std::vector < std::vector <double> > &particles);
  void ListFiles(void) const;
  std::string GetCurrentFileName(void) const { return fcname; }
  int SetFileType(std::string ftype);

protected:  
  int ReadRecord(std::vector < std::vector <double> > &particles);  
  int ReadRecordOut(std::vector < std::vector <double> > &particles);

protected:  
  std::vector<std::string> fname_list;
  std::string              fcname;
  std::fstream             fdata;
  int                      fid;
  int                      n_columns, skip_lines;
  int                     (FileData::*freadfn)(std::vector < std::vector <double> > &particles);
  int                      debugl;
};


#endif    // 

