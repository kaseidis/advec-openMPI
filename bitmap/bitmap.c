#include "bitmap.h"
#include <stdlib.h>

// Adapted from https://gist.github.com/guohai/5848088
int save_bitmap(RgbTriple* data, int32_t width, int32_t height, const char* filename)
{
  uint32_t data_size = width * height * sizeof(RgbTriple);

  BitmapFileHeader fhead;
  fhead.type = 0x4D42; // 'BM'
  fhead.size = sizeof(BitmapFileHeader) + sizeof(BitmapInfoHeader) + data_size;
  fhead.reserved1 = 0;
  fhead.reserved2 = 0;
  fhead.off_bits = sizeof(BitmapFileHeader) + sizeof(BitmapInfoHeader);

  BitmapInfoHeader finfo;
  finfo.size = sizeof(BitmapInfoHeader);
  finfo.width = width;
  finfo.height = height;
  finfo.planes = 1;
  finfo.bit_count = sizeof(RgbTriple) * 8;
  finfo.compression = 0;
  finfo.size_image = data_size;
  finfo.x_pels_per_meter = 0;
  finfo.y_pels_per_meter = 0;
  finfo.clr_used = 0;
  finfo.clr_important = 0;

  FILE* fp = fopen(filename, "wb");
  if (fp) {
    fwrite(&fhead, 1, sizeof(BitmapFileHeader), fp);
    fwrite(&finfo, 1, sizeof(BitmapInfoHeader), fp);
    fwrite(data, width * height * sizeof(RgbTriple), 1, fp);
    fclose(fp);
    return 0;
  } else {
    fprintf(stderr, "ERROR: Couldn't open bitmap file %s for writing\n", filename);
    return 1;
  }
}

RgbTriple* load_bitmap(const char* filename,
                       BitmapFileHeader* fhead,
                       BitmapInfoHeader* finfo)
{
  FILE* f = fopen(filename, "rb");
  if (f == NULL) {
    fprintf(stderr, "ERROR: Couldn't read bitmap file: %s\n", filename);
    return NULL;
  }
  fread(fhead, sizeof(BitmapFileHeader), 1, f);
  fread(finfo, sizeof(BitmapInfoHeader), 1, f);
  if (finfo->bit_count != sizeof(RgbTriple) * 8) {
    fprintf(stderr,
            "ERROR: Attempted to read RgbTriple data from %s but file uses wrong pixel "
            "bit count (%d bits)\n",
            filename,
            (finfo->bit_count));
    return NULL;
  }
  fseek(f, fhead->off_bits, SEEK_SET);
  RgbTriple* img = malloc(finfo->size_image);
  fread(img, finfo->size_image, 1, f);
  return img;
}
