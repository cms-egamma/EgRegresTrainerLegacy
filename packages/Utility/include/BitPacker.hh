
#include <stdint.h>
class BitPacker {
  BitPacker()=delete;
  ~BitPacker()=delete;
public:
  
  template<typename T>
  static T unpack(const T packedValue, const uint32_t offset, const uint32_t width){
    return (packedValue >> offset) & ((1 << width) - 1);
  }
  
  template<typename T>
  static T pack(T packedValue,const T value, const uint32_t offset, const uint32_t width) {
    const T mask = ((1 << width) - 1) << offset;
    packedValue &= ~mask;
    packedValue |= (value & ((1U << width) - 1)) << offset;
    return packedValue;
  }

  

};
