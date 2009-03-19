/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkTipsyReader.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkTipsyReader - Read ASCII or binary particle 
//                            data and (optionally) one scalar 
//                            value associated with each particle.
// .SECTION Description
// vtkTipsyReader reads either a binary or a text file of 
//  particles. Each particle can have associated with it an optional
//  scalar value. So the format is: x, y, z, scalar 
//  (all floats or doubles). The text file can consist of a comma 
//  delimited set of values. In most cases vtkTipsyReader can 
//  automatically determine whether the file is text or binary. 
//  The data can be either float or double. 
//  Progress updates are provided. 
//  With respect to binary files, random access into the file to read 
//  pieces is supported.
//  

#ifndef __vtkTipsyReader_h
#define __vtkTipsyReader_h

#include "vtkPolyDataAlgorithm.h"

#define VTK_FILE_BYTE_ORDER_BIG_ENDIAN 0
#define VTK_FILE_BYTE_ORDER_LITTLE_ENDIAN 1


class VTK_IO_EXPORT vtkTipsyReader : public vtkPolyDataAlgorithm
{
public:
  static vtkTipsyReader *New();
  vtkTypeRevisionMacro(vtkTipsyReader,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);   

  // Description:
  // Specify file name.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // These methods should be used instead of the SwapBytes methods.
  // They indicate the byte ordering of the file you are trying
  // to read in. These methods will then either swap or not swap
  // the bytes depending on the byte ordering of the machine it is
  // being run on. For example, reading in a BigEndian file on a
  // BigEndian machine will result in no swapping. Trying to read
  // the same file on a LittleEndian machine will result in swapping.
  // As a quick note most UNIX machines are BigEndian while PC's
  // and VAX tend to be LittleEndian. So if the file you are reading
  // in was generated on a VAX or PC, SetDataByteOrderToLittleEndian 
  // otherwise SetDataByteOrderToBigEndian. Not used when reading
  // text files. 
  void SetDataByteOrderToBigEndian();
  void SetDataByteOrderToLittleEndian();
  int GetDataByteOrder();
  void SetDataByteOrder(int);
  const char *GetDataByteOrderAsString();

  // Description:
  // Set/Get the byte swapping to explicitly swap the bytes of a file.
  // Not used when reading text files.
  vtkSetMacro(SwapBytes,int);
  int GetSwapBytes() {return this->SwapBytes;}
  vtkBooleanMacro(SwapBytes,int);
  
  // Description:
  // Default: 1. If 1 then each particle has a value associated with it.
  vtkSetMacro(HasScalar,int);
  vtkGetMacro(HasScalar,int);
  vtkBooleanMacro(HasScalar,int);

  // Description:
  // Get/Set the file type.  The options are:
  // - FILE_TYPE_IS_UNKNOWN (default) the class 
  //     will attempt to determine the file type.
  //     If this fails then you should set the file type
  //     yourself.
  // - FILE_TYPE_IS_TEXT the file type is text.
  // - FILE_TYPE_IS_BINARY the file type is binary.
  vtkSetClampMacro(FileType, int, FILE_TYPE_IS_UNKNOWN, FILE_TYPE_IS_TIPSY_NATIVE);
  vtkGetMacro(FileType, int);
  void SetFileTypeToUnknown() {this->SetFileType(FILE_TYPE_IS_UNKNOWN);}
  void SetFileTypeToStandard() {this->SetFileType(FILE_TYPE_IS_TIPSY_STANDARD);}
  void SetFileTypeToNative() {this->SetFileType(FILE_TYPE_IS_TIPSY_NATIVE);}

#if 0
  // Description:
  // Get/Set the data type.  The options are:
  // - VTK_FLOAT (default) single precision floating point.
  // - VTK_DOUBLE double precision floating point.
  vtkSetClampMacro(DataType, int, VTK_FLOAT, VTK_DOUBLE);
  vtkGetMacro(DataType, int);
  void SetDataTypeToFloat() {this->SetDataType(VTK_FLOAT);}
  void SetDataTypeToDouble() {this->SetDataType(VTK_DOUBLE);}
#endif


protected:
  vtkTipsyReader();
  ~vtkTipsyReader();

  void OpenFile();

  char *FileName;
  ifstream *File;

  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  
  // Description:
  // Determine the type of file based on an analysis of its contents.
  // Up to 5000 bytes of the file are read and classified. The classification
  // of a file as either binary or text is based on the proportions of bytes in
  // various classifications. The classification of the file is not infallible
  // but should work correctly most of the time. If it fails, use SetFileTypeToText()
  // or SetFileTypeToBinary() to set the file type.
  // This algorithm probably only identifies ASCII text correctly and will not 
  // work for UTF-8 UCS-2 (or UTF-16) or UCS-4 or EBCIDIC.
  int DetermineFileType();
  
  // Description:
  // Update of the progress.
  void DoProgressUpdate( size_t & bytesRead, size_t & fileLength );

  //BTX
  // Description:
  // Enumerate the supported file types.
  // <pre>
  // - FILE_TYPE_IS_UNKNOWN, (default) the class will attempt to determine the file type.
  // - FILE_TYPE_IS_TIPSY_STANDARD, the file type is Tipsy standard.
  // - FILE_TYPE_IS_TIPSY_NATIVE, the file type is Tipsy native.
  // </pre>
  enum FILE_TYPE { FILE_TYPE_IS_UNKNOWN = 0, 
    FILE_TYPE_IS_TIPSY_STANDARD, FILE_TYPE_IS_TISPY_NATIVE };
  //ETX
  
  int HasScalar;
  // Description:
  // Used to decide which reader should be used.
  int FileType;
  // Description:
  // Used to specify the data type.
  int DataType;

  // Description:
  // Set an alliquot of bytes.  
  size_t Alliquot;
  // Description:
  // Count of the number of alliquots processed.
  size_t Count;

  int SwapBytes;
  size_t NumberOfPoints;

private:
  vtkTipsyReader(const vtkTipsyReader&);  // Not implemented.
  void operator=(const vtkTipsyReader&);  // Not implemented.
};

#endif
