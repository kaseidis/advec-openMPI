#ifndef MPCS_HPC_WINTER_2021_SCRATCH_BITMAP_BITMAP_H
#define MPCS_HPC_WINTER_2021_SCRATCH_BITMAP_BITMAP_H

#include <stdint.h>
#include <stdio.h>

// References:
// * https://docs.microsoft.com/en-us/windows/win32/gdi/bitmap-storage
// * https://docs.microsoft.com/en-us/windows/win32/gdi/bitmap-structures
// * https://docs.microsoft.com/en-us/windows/win32/winprog/windows-data-types
// * https://gist.github.com/guohai/5848088

// See:
// https://docs.microsoft.com/en-us/windows/win32/api/wingdi/ns-wingdi-bitmapfileheader
typedef struct {
  uint16_t type;      // The file type; must be BM
  uint32_t size;      // The size, in bytes, of the bitmap file.
  uint16_t reserved1; // Reserved; must be zero.
  uint16_t reserved2; // Reserved; must be zero.
  // The offset, in bytes, from start of BITMAPFILEHEADER structure to bitmap bits.
  uint32_t off_bits;
} __attribute__((packed)) BitmapFileHeader;

// See: https://docs.microsoft.com/en-us/previous-versions//dd183376(v=vs.85)
typedef struct {
  // The number of bytes required by the structure.
  uint32_t size;
  // The width of the bitmap, in pixels.
  int32_t width;
  // The height of the bitmap, in pixels.
  // * If biHeight is positive, its origin is the lower-left corner.
  // * If biHeight is negative, its origin is the upper-left corner
  int32_t height;
  // The number of planes for the target device. Must be 1
  uint16_t planes;
  // The number of bits-per-pixel.
  uint16_t bit_count;
  // Type of compression.  "BI_RGB" for uncompressed
  uint32_t compression;
  // The size, in bytes, of the image. Can be 0 for BI_RGB (uncompressed) bitmaps
  uint32_t size_image;
  // The horizontal resolution, in pixels-per-meter, of the target device. Can be 0
  int32_t x_pels_per_meter;
  // The vertical resolution, in pixels-per-meter, of the target device.
  int32_t y_pels_per_meter;
  // The number of color indices used by the bitmap. If 0, inferred from biBitCount
  uint32_t clr_used;
  // The number of color indexes for displaying the bitmap.  If 0, inferred
  uint32_t clr_important;
} __attribute__((packed)) BitmapInfoHeader;

// See: https://docs.microsoft.com/en-us/windows/win32/api/wingdi/ns-wingdi-rgbtriple
typedef struct {
  unsigned char blue;
  unsigned char green;
  unsigned char red;
} RgbTriple;

int save_bitmap(RgbTriple* data, int32_t width, int32_t height, const char* filename);
RgbTriple* load_bitmap(const char* filename,
                       BitmapFileHeader* fhead,
                       BitmapInfoHeader* finfo);

#endif // MPCS_HPC_WINTER_2021_SCRATCH_BITMAP_BITMAP_H
