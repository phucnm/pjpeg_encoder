/* input_stream.hpp
   CSC 485B/CSC 578B/SENG 480B - Summer 2020

   Definition of a bitstream class which complies with the bit ordering
   required by the gzip format.

   (The member function definitions are all inline in this header file 
    for convenience, even though the use of such long inlined functions
    might be frowned upon under some style manuals)

   B. Bird - 06/19/2020
*/ 

#ifndef INPUT_STREAM_HPP
#define INPUT_STREAM_HPP

#include <iostream>
#include <cstdint>

/* These definitions are more reliable for fixed width types than using "int" and assuming its width */
using u8 = std::uint8_t;
using u16 = std::uint16_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;



class InputBitStream{
public:
    /* Constructor */
    InputBitStream( std::istream& input_stream ): bitvec{0}, numbits{8}, infile{input_stream}, done{false}, last_real_bit{0} {

    }

    /* Destructor */
    virtual ~InputBitStream(){
    }

    /* Read an entire byte from the stream, with the least significant bit read first */
    unsigned char read_byte(){
        return read_bits(8);
    }

    /* Read a 32 bit unsigned integer value (LSB first) */
    u32 read_u32(){
        return read_bits(8) | (read_bits(8)<<8) | (read_bits(8)<<16) | (read_bits(8)<<24);
    }

    /* Read a 16 bit unsigned short value (LSB first) */
    u16 read_u16(){
        return read_bits(8) | (read_bits(8)<<8);
    }

    /* Read the lowest order num_bits bits from the stream into a u32,
       with the least significant bit read first.
    */
    u32 read_bits(int num_bits){
        u32 result{};
        for (int i = 0; i < num_bits; i++)
            result |= read_bit()<<i;
        return result;
    }

    /* Read a single bit b (stored as the LSB of an unsigned int)
       from the stream */
    unsigned int read_bit(){
        //This has been set up to emit an infinite number of copies
        //of the last bit once EOF is reached (so if the last bit
        //in the file is a 1, any subsequent call to read_bit will
        //return 1).
        if (numbits == 8)
            input_byte();
        if (!done)
            last_real_bit = (bitvec>>(numbits++))&0x1;
        
        return last_real_bit;
    }

    /* Flush the currently stored bits*/
    void flush_to_byte(){
        numbits = 8; //Force the next read to read a byte from the input file
    }
private:
    void input_byte(){
        char c;
        if (!infile.get(c)){
            done = true;
            numbits = 0;
            return;
        }
        bitvec = (unsigned char)c;
        numbits = 0;
    }
    u32 bitvec;
    u32 numbits;
    std::istream& infile;
    bool done;
    unsigned int last_real_bit;
};


#endif 