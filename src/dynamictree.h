#ifndef _DYNAMIC_TREE_H
#define _DYNAMIC_TREE_H

#include <string>

#include <TMath.h>

/*
 * SOME CONSTANTS
 */

#define EVENT_HEADER_HEADER 0x08  // 0b1000
#define EVENT_UUID_HEADER   0x08  // 0b1000
#define WIRE_ENTRY_HEADER   0X00  // 0b0000
#define SCINT_ENTRY_HEADER  0x07  // 0b0111
#define END_OF_FILE_HEADER  0x0F  // 0b1111

#define MASK_HEADER_SHIFT 28
#define MASK_NSCINT_SHIFT 24
#define MASK_NWIRES_SHIFT 16
#define MASK_DC_CHN_SHIFT 16
#define MASK_SC_CHN_SHIFT 16
#define MASK_BLK_ID_SHIFT 16
#define MASK_BLK_CK_SHIFT 24

#define MASK_HEADER      0xF0000000
#define MASK_NSCINT      0x0F000000
#define MASK_NWIRES      0x00FF0000
#define MASK_PATTRN      0x000000FF
#define MASK_EVENT_UUID  0x0000FFFF
#define MASK_BLOCK_UUID  0x00FF0000
#define MASK_DC_CHN      0x007F0000
#define MASK_DC_TIME     0x0000FFFF
#define MASK_SCINT_CHN   0x00070000
#define MASK_SCINT_TIME  0x0000FFFF
#define MASK_TOT_EVENTS  0x00FFFFFF
#define MASK_BLOCK_CHK   0x01000000

/*
 * USEFUL C MACROS 
 */

#define GET_TOKEN(data, mask, nbits) (data & mask) >> nbits

#define GET_HEADER(X) GET_TOKEN(X, MASK_HEADER, MASK_HEADER_SHIFT)
#define GET_NSCINT(X) GET_TOKEN(X, MASK_NSCINT, MASK_NSCINT_SHIFT)
#define GET_NWIRES(X) GET_TOKEN(X, MASK_NWIRES, MASK_NWIRES_SHIFT)
#define GET_PATTRN(X) (X & MASK_NWIRES)

#define GET_DC_CHANNEL(X) GET_TOKEN(X, MASK_DC_CHN, MASK_DC_CHN_SHIFT)
#define GET_DC_TIME(X) (X & MASK_DC_TIME)

#define GET_SCINT_CHANNEL(X) GET_TOKEN(X, MASK_SCINT_CHN, MASK_SC_CHN_SHIFT)
#define GET_SCINT_TIME(X) (X & MASK_DC_TIME)

#define GET_EVENT_UUID(X) (X & MASK_EVENT_UUID)
#define GET_BLOCK_UUID(X) GET_TOKEN(X, MASK_BLOCK_UUID, MASK_BLK_ID_SHIFT)

#define GET_BLOCK_CHK(X) GET_TOKEN(X, MASK_BLOCK_CHK, MASK_BLK_CK_SHIFT)

#define GET_TOTAL_EVENTS(X) (X & MASK_TOT_EVENTS)

//#define CLR_HEAP_ARR(X) std::fill(X, X+sizeof(X), 0)
#define CLR_HEAP_ARR(X) memset(X, 0, sizeof(X))

/*
 * MESSAGES DECLARATIONS
 */

extern const char* E_CORR_ENTRY;
extern const char* E_FILE_OPEN;
extern const char* E_MEM_ERROR;
extern const char* E_WRNG_ENTRY;
extern const char* I_EVENT_INFO;
/*
 * FUNCTION DECLARATIONS
 */

void parse(std::string in_file_url);
void parseAll();
void hist();
Double_t mygaus(Double_t *x, Double_t *par);

#endif
