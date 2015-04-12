//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// BitChar - Fixed-size bitset of 8 bits
//

#include <iostream>
#include <assert.h>
#include "BitChar.h"

static unsigned char bc_mask[]={0x80,0x40,0x20,0x10,0x08,0x04,0x02,0x01};

//
void BitChar::set(unsigned char idx, bool v)
{
    if(v)
        m_data = m_data | bc_mask[idx];
    else
        m_data = m_data & ~bc_mask[idx];
}

//
bool BitChar::test(unsigned char idx) const
{
    return m_data & bc_mask[idx];
}

//
void BitChar::flip(unsigned idx)
{
    m_data = m_data ^ bc_mask[idx];
}

//
void printBinary(std::ostream& out, const BitChar& bc)
{
    for(int i = 7; i >= 0; i--)
    {
        if(bc.test(i))
            out << "1";
        else
            out << "0";
    }
}

//
std::ostream& operator<<(std::ostream& out, const BitChar& bc)
{
    int temp = bc.m_data;
    out << temp;
    return out;
}

//
std::istream& operator>>(std::istream& in, BitChar& bc)
{
    int temp;
    in >> temp;
    bc.m_data = (unsigned char)temp;
    return in;
}

// 
void BitChar::write(std::ostream& out)
{
    out.write((char*)&m_data, sizeof(m_data));
}

//
void BitChar::read(std::istream& in)
{
    in.read((char*)&m_data, sizeof(m_data));
}
