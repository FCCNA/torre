#include <iostream>
#include <fstream>

#include <TTree.h>
#include <TFile.h>

#include "dynamictree.h"

/*
 * MESSAGES DEFINITIONS
 */

const char* E_CORR_ENTRY = "Big Jupiter! There is a %s outside an event!\n"
                           "This usually just means that the input file is "
                           "corrupted and the previous event is garbage!";

const char* E_FILE_OPEN = "Err... I cannot open the input file %s, are you "
                          "sure that it *actually* exists?";

const char* E_MEM_ERROR = "Cannot create %s: not enough memory!";

const char* E_WRNG_ENTRY = "\nI expected this was a %s, but it seems I was "
                           "wrong...\nActually it's more probable the input "
                           "file is just corrupted.\nI'm bailing out, you are"
                           " out of luck, let's try again after you fixed the"
                           " input file ;)\n";

const char* I_EVENT_INFO = "This was block %d event %d [%ld, %ld]";

/*
 * FUNCTIONS DEFINITIONS
 */

using namespace std;

int main()
{
  /*
	*
	* Here list the name of the file to be parsed
	*/

	parse("data/Tracker_231122_1212.txt");
  return 0;
}

void parse(std::string in_file_url)
{
  /*
   * ROOT macro for parsing plain-text encoded data files
   *
   * @desc: This macro reads the user provided input file, then creates
   *        a TTree object from the its content and finally saves the TTree
   *        to a file named "data_tree.root". If the destination file already
   *        exists, it will be deleted and recreated with the new TTree object.
   *        The input file format is described in the comments below.
   */

  string out_file_name = in_file_url + string(".root");

  unsigned long num_events_read = 0; // this variable holds the number of
                                     // evets that has been actually read
                                     // from the input file.

  bool end_flag = false; // flag used to check if the input file is good

  unsigned int num_of_scint = 0; // this variable holds the nubmer of
                                 // scintillators which have been triggered
  unsigned int num_of_wires = 0; // this variable holds the number of
                                 // drift chambers which have been triggered

  float scint_time[8];  // buffer used to store the times of the
                        // scintillators.
  float wire_time[96]; // buffer used to store the times of the
                        // drift chambers.

  int scint_index[8];  // array used to store the indexes of the
                       // scintillators that have been triggered
  int wire_index[96]; // array used to store the indexes of the
                       // driftchambers that have been triggered

  unsigned int event_uuid = 0;  // variable to hold the current event UUID
  unsigned int block_uuid = 0;  // variable to hold the current block UUID
  unsigned int block_count = 0; // variable to hold the numbers of blocks
                                // of events that has been read.

  std::string line1; // string used to store the first line of the input file
  std::string line2; // string used to store the sencond line of the input file

  std::ifstream in_file; // the input file itself

  in_file.open(in_file_url);
  // check if the input file is actually open...
  if (!in_file.is_open())
  {
    // ...and if it isn't, just print a pretty error message.
    Error("", E_FILE_OPEN, in_file_url.c_str());
    return;
  }

  // If everything has gone as expected, now we will create a TTree object
  // and also a TFile object to save the former to a file.
  TFile out_file(out_file_name.c_str(), "recreate");
  TTree *tdata = new TTree("tdata", "data");

  // Again, if there is enough free memory to create the TTtree...
  if (tdata == NULL)
  {
    // ...close everything and gently inform the user.
    Error("", E_MEM_ERROR, "the main TTree [tdata]");
    in_file.close();
    out_file.Close();
    return;
  }

  // creating the branches of the main TTree and set to each its on
  // associated variable
  tdata->Branch("SCINT_COUNT", &num_of_scint, "SCINTS_COUNT/I");
  tdata->Branch("SCINT_INDEX", &scint_index, "SCINT_INDEX[SCINTS_COUNT]/I");
  tdata->Branch("SCINT_TIME", &scint_time, "SCINT_TIME[SCINTS_COUNT]/F");

  tdata->Branch("WIRES_COUNT", &num_of_wires, "WIRES_COUNT/I");
  tdata->Branch("WIRE_INDEX", &wire_index, "WIRE_INDEX[WIRES_COUNT]/I");
  tdata->Branch("WIRE_TIME", &wire_time, "WIRE_TIME[WIRES_COUNT]/F");

  /*
   *  NOTE: A BRIEF DESCRIPTION OF THE INPUT FILE STRUCTURE
   *
   *  1. The basic file structure
   *
   *  The input file is basically a plain-text sqeuence of entries. Each entry
   *  is a 32bit binary encoded word which contains different informations
   *  depending of the type of the entry itself.
   *
   *  The following example explains better how the input file is organized: 
   *
   *  +-----------------------------------------+
   *  |  Wednesday, December 02, 2015 11:37 AM  | <--[FILE_HEADER]
   *  |  Tracking system                        |
   *  +-----------------------------------------+                          ---
   *  |  20151202 1137                          | <--[BLOCK_HEADER]         |
   *  +-----------------------------------------+                  ---      |
   *  |  -2147483648                            | <--[ EVENT_UUID ] |       |
   *  |  -2113470392                            | <--[EVENT_HEADER] |       |
   *  +-----------------------------------------+                   |       |
   *  |  918832                                 |                   |       |
   *  |  2229615                                |                   E       |
   *  |  2950741                                |                   V       |
   *  |  3212189                                | <--[   WIRES    ] E       |
   *  |  3933655                                |    [  ENTRIES   ] N       |
   *  |  5178074                                |                   T       B
   *  |  5899934                                |                   |       L
   *  +-----------------------------------------+                   |       O
   *  |  1879246513                             | <--[   SCINT    ] |       C
   *  |  1879442461                             |    [  ENTRIES   ] |       K
   *  +-----------------------------------------+                  ---
   *     ....  OTHER EVENTS ...                                             O
   *  +-----------------------------------------+                  ---      F
   *  | -2147482647                             | <--[ EVENT UUID ] |
   *  | -2113404856                             | <--[EVENT HEADER] |       E
   *  +-----------------------------------------+                   |       V
   *  |  263732                                 |                   |       E
   *  |  983765                                 |                   E       N
   *  |  2229809                                |                   V       T
   *  |  2949776                                | <--[   WIRES    ] E       S
   *  |  3539992                                |    [  ENTRIES   ] N       |
   *  |  4261123                                |                   T       |
   *  |  5244261                                |                   |       |
   *  |  5964813                                |                   |       |
   *  +-----------------------------------------+                   |       |
   *  |  1879246513                             | <--[   SCINT    ] |       |
   *  |  1879442077                             |    [  ENTRIES   ] |       |
   *  +-----------------------------------------+                  ---     ---
   *    ... OTHER BLOCKS ...       
   *  +-----------------------------------------+
   *  |  -16767216                              | <--[END_OF_FILE]
   *  +-----------------------------------------+
   *
   *  The first two lines of the files contain the date and time the file
   *  itself was created and the software used to collect the data.
   *
   *  Then several BLOCKS OF EVENTS follows. Each BLOCK starts with a
   *  BLOCK_HEADER that contains the date and time in which the first
   *  event of the block, then there are several EVENTs.
   *
   *  Each EVENT section starts with the EVENT_UUID, which is a monotonic
   *  increasing number that identifies the event. Then there is the
   *  EVENT_HEADER which contains information about the number of
   *  scintillators and drift chambers that have detected the event.
   *  Then there are as many entries as number of drift chambers that
   *  detected the event, and finally as many entries as number of
   *  scintillators that detected the event.
   *
   *  Finally there is the END_OF_FILE entry which contains the total number
   *  ov events that have been recorderd in the file.
   *
   *  2. Entries structure
   *
   *    As previously said, each entry is basically a 32bit word. The last
   *    four bits are the ENTRY_HEADER that identifies the type of the entry,
   *    as you can see from the following table:
   *
   *        +--------------------------------+
   *        |          ENTRY HEADERS         |
   *        +--------+------+----------------+
   *        | BINARY |  HEX |   ENTRY TYPE   |
   *        +--------+------+----------------+
   *        |  1000  | 0x08 |  EVENT_HEADER  |
   *        +--------+------+----------------+
   *        |  1000  | 0x08 |    EVENT_UUID  |
   *        +--------+------+----------------+
   *        |  0000  | 0x00 |    WIRE_ENTRY  |
   *        +--------+------+----------------+
   *        |  0111  | 0x07 |   SCINT_ENTRY  |
   *        +--------+------+----------------+
   *        |  1111  | 0x0F |   END_OF_FILE  |
   *        +--------+------+----------------+
   *
   *    2.3 EVENT_HEADER
   *
   *      |32  29|28          25|24          17|16       9|8       1|
   *      +------+--------------+--------------+----------+---------+
   *      | 1000 | num_of_scint | num_of_wires | 00000000 | pattern |
   *      +------+--------------+--------------+----------+---------+
   *
   *      - bits   1-8: are the scintillator pattern: each bit corresponds to
   *                    a scintillator (I need more information, the manual is
                        not useful).
   *      - bits  9-16: are unused and should be all 0.
   *      - bits 17-24: encode the nunmber of drift chambers (num_of_wires)
   *                    that detected the event.
   *      - bits 25-28: encode the number of scintillators (num_of_scint)
   *                    that detected the particle.
   *      - bits 29-32: are the ENTRY_HEADER.
   * 
   *    2.4 EVENT_UUID
   *
   *      NOTE: This entry is missing from the file which describes how
   *            the encoding is done, so what follows has been deduced
   *            by reverse-engineering and thus it may be incorrect...
   *            you have been warned ;)
   *
   *      |32  29|28                  17|16                        1|
   *      +------+----------------------+---------------------------+
   *      | 1000 |      block_uuid      |        event_uuid         |
   *      +------+----------------------+---------------------------+
   *
   *      - bits  1-16: encode the monotonic event number for the
   é                    current event block (event_uuid).
   *      - bits 17-28: encode the monotonic block number (block_uuid).
   *      - bits 29-32: are the ENTRY_HEADER.
   *
   *    2.5 WIRE_ENTRY
   *
   *      |32  29|28   25|  24 |23          17|16                  1|
   *      +------+-------+-----+--------------+---------------------+
   *      | 0000 |       |  0  | wire_channel |      wire_time      |
   *      +------+-------+-----+--------------+---------------------+
   *
   *      - bits  1-16: encode the time at which the drift chamber detected
   *                    the particle, in units of 1.04ns (wire_time).
   *      - bits 17-23: encodes the drift chamber index (wire_channel).
   *      - bit     24: is unused and should be 0.
   *      - bits 25-28: are unused.
   *      - bits 29-32: are the ENTRY_HEADER.
   *
   *    2.6 SCINT_ENTRY
   *
   *      |32  29|28   21|  20 |19          17|16                  1|
   *      +------+-------+-----+--------------+---------------------+
   *      | 0111 |       |  0  |scint_channel |      scint_time     |
   *      +------+-------+-----+--------------+---------------------+
   *
   *      - bits  1-16: encode the time at which the scintillator detected
   *                    the particle, in unit of 1.04ns (scint_time).
   *      - bits 17-19: encodes the scintillator index (scint_channel).
   *      - bit     20: is unused and should be 0.
   *      - bits 21-28: are unused.
   *      - bits 29-32: are the ENTRY_HEADER.
   *
   *    2.6 END_OF_FILE
   *
   *      |32  29|28   25|24                                       1|
   *      +------+-------+------------------------------------------+
   *      | 1111 |  1111 |               num_of_events              |
   *      +------+-------+------------------------------------------+
   *
   *      - bits  1-24: encode the number of events (num_of_events)
   *      - bits 25-28: are unused and must be all 1.
   *      - bits 29-32: are the ENTRY_HEADER.
   *
   *  For a full description of the file structure, please read the file
   *  'DAQ for Tracking system.pdf'
   */


  // reading the FILE_HEADER...
  getline(in_file, line1);
  getline(in_file, line2);

  // ... and printing it to the console
  std::cout << std::endl;
  std::cout << "Acquisition date: " << line1 << std::endl;
  std::cout << "Software: " << line2 << std::endl << std::endl;

  // Very well, it's time to read the input file and populate the TTree.

  while(!end_flag && in_file.good())
  {
    unsigned long entry;
    unsigned int entry_header;

    string entry_line;

    in_file >> entry; // reading the next 
    entry_header = GET_HEADER(entry);

    /*
     * Let's try to explain what we are doing: we already read the
     * FILE_HEADER so from now on there will be only BLOCK_HEADERs,
     * EVENT_HEADERs or one of the other entries. However, the EVENT_HEADERs
     * already contain the number of SCINT entries and WIRE entries for that
     * event and thus we will use two for loop to read them.
     *
     * Now, since all the SCINT and WIRE entries have been read in the for
     * loops, here we can read a BLOCK_HEADER, an EVENT_UUID, or an
     * END_OF_FILE entry.
     *
     * Of course, this is not true if the file is somehow corrupted
     * (for example the EVENT_HEADER says that only two scintillators
     * have been triggered in the current event, but there are actually
     * three SCINT entries... and things like that, unfortunately, happens!)
     * In such cases, the only thing to do is to blame LabView and/or who 
     * created the data acquisition program...
     */

    switch(entry_header)
    {
      case EVENT_UUID_HEADER:
      {
        // this is the most common case, we have found and EVENT, which
        // always starts with and EVENT_UUID entry. Next there is the
        // event header. 

        unsigned int header;
        // unsigned char pattern;  // unused
        unsigned long event_header;

        event_uuid = GET_EVENT_UUID(entry);
        block_uuid = GET_BLOCK_UUID(entry);

        // read the EVENT_HEADER entry
        in_file >> event_header;

        // Uncomment the following lines when you need to debug the code.
        //Info("", "Reading block %d event %d", block_uuid, event_uuid);

        // read the event properties
        header       = GET_HEADER(event_header);
        num_of_scint = GET_NSCINT(event_header);
        num_of_wires = GET_NWIRES(event_header);
        // pattern      = GET_PATTRN(event_header);  // unused

        // Uncomment the following lines if you want a very
        // very very very very verbose output.
        //Info("", "Found new event: %u scints %u wires",
        //         num_of_scint, num_of_wires);

        if(header != EVENT_HEADER_HEADER)
        {
          Error("", E_WRNG_ENTRY, "EVENT_HEADER");
          end_flag=true;
          break;
        }

        if ((num_of_scint >= 8) && (num_of_wires >= 96))
        {
          end_flag=true;
          Warning("", "Wrong number of wires/scint!");
        }

        // Clean the arrays
        CLR_HEAP_ARR(scint_index);
        CLR_HEAP_ARR(scint_time);
        CLR_HEAP_ARR(wire_index);
        CLR_HEAP_ARR(wire_time);

        // read all the WIRE entries for the current event
        for(unsigned int i = 0; (!end_flag) && (i < num_of_wires); i++)
        {
          unsigned long entry_wire;
          unsigned int wire_header;

          // read the WIRE_ENTRY
          in_file >> entry_wire; 

          // retrieve the ENTRY_HEADER
          wire_header = GET_HEADER(entry_wire);

          // check the type of the ENTRY is correct...
          if(wire_header != WIRE_ENTRY_HEADER)
          {
            // ...otherwise print an error message
            Error("", E_WRNG_ENTRY, "WIRE_ENTRY");
            Info("", I_EVENT_INFO,
                     block_uuid, event_uuid,
                     entry, event_header);
            end_flag=true;
            break;
          }

          // store the entry properties
          wire_index[i] = GET_DC_CHANNEL(entry_wire);
          wire_time[i] = GET_DC_TIME(entry_wire);
        } 

        // read all the SCINT entries for the current event
        for(unsigned int i = 0; (!end_flag) && (i < num_of_scint); i++)
        {
          unsigned long entry_scint;
          unsigned int scint_header;

          // read the SCINT_ENTRY
          in_file >> entry_scint;

          // retrieve the ENTRY_HEADER
          scint_header = GET_HEADER(entry_scint);

          // check the type of the ENTRY is correct...
          if(scint_header != SCINT_ENTRY_HEADER)
          {
            // ...otherwise print an error message
            Error("", E_WRNG_ENTRY, "SCINT_ENTRY");
            Info("", I_EVENT_INFO, 
                     block_uuid, event_uuid,
                     entry, event_header);
            end_flag=true;
            break;
          }

          // store the entry properties
          scint_index[i] = GET_SCINT_CHANNEL(entry_scint);
          scint_time[i] = GET_SCINT_TIME(entry_scint);
        }

        // If everything is ok...
        if (!end_flag)
        {
          // ...fill the TTree and incremente the number of events read.
          tdata->Fill();
          num_events_read++;
        }
        break;
      }
      case END_OF_FILE_HEADER:
      {
        // we have reached the end of the file...
        // let's do an overall sanity check!
        unsigned long tot_events = GET_TOTAL_EVENTS(entry);
        float corr_perc = 100.0 * (1.0 - (num_events_read/tot_events));

        Info("", "Reached the end of the file %s\n",
                 in_file_url.c_str());
        Info("", "Total block:\t%8u", block_count);
        Info("", "Total events expected:\t%8ld", tot_events);
        Info("", "Total events read:\t%8ld", num_events_read);
        Info("", "Data corruption is below %d%%\n",
                  (unsigned int) std::ceil(corr_perc));
        // Oh! And remember to set the end_flag!
        end_flag = true;
        break;
      }
      case SCINT_ENTRY_HEADER:
      {
        // This will be execute only if the input file is corrupted:
        // more precisely when the EVENT_HEADER of previous even lied
        // about the number of SCINT entries.
        // 
        // Shame on LabView (see the NOTE above)!
        Warning("", E_CORR_ENTRY, "SCINT entry");
        Info("", I_EVENT_INFO, block_uuid, event_uuid, 0, 0);

        break;
      }
      case WIRE_ENTRY_HEADER:
      {
        // This will be execute only if the input file is corrupted:
        // more precisely when the EVENT_HEADER of the previous event
        // lied about the number of WIRE entries.
        // 
        // Shame on LabView (see the NOTE above)!

        // NOTE: the problem here is that a WIRE_ENTRY and a BLOCK_HEADER
        //       has the same ENTRY_HEADER... However, the bit 24 is zero
        //       for a WIRE_ENTRY but it is not forr a BLOCK_HEADER.
        //       Banzai!
        //
        //       This check is preformed in the 'default' case.
      }
      default:
      {
        if (GET_BLOCK_CHK(entry) == 0)
        {
          Warning("", E_CORR_ENTRY, "WIRE entry");
          Info("", I_EVENT_INFO, block_uuid, event_uuid, 0, 0);
        }
        else
        {
          // this entry can only be a BLOCK_HEADER one. We have already read
          // the date, now there should be the time.
          unsigned long time;
          in_file >> time;
          Info("", "Found new block of events [%ld - %ld]", entry, time);
          block_count++;
          break;
        }
      }
    }
  }

  // That's all, now we just have to save our TTree to a file.
  in_file.close();
  tdata->Write();
  out_file.Close();

  return;
}
