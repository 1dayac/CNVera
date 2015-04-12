//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// MultiOverlap.h - Data structure containing a set
// of overlaps for a given read
//
#include <algorithm>
#include <iostream>
#include "MultiOverlap.h"
#include "Alphabet.h"

MultiOverlap::MultiOverlap(const std::string& rootID, 
                           const std::string& rootSeq,
                           const std::string rootQual) : m_rootID(rootID), 
                                                         m_rootSeq(rootSeq),
                                                         m_rootQual(rootQual)
{
    
}

//
void MultiOverlap::add(const std::string& seq, const Overlap& ovr)
{
    MOData mod(seq, ovr);

    // Swap root read into first position if necessary
    if(mod.ovr.id[0] != m_rootID)
        mod.ovr.swap();
    assert(mod.ovr.id[0] == m_rootID);

    // RC the sequence if it is different orientation than the root
    if(mod.ovr.match.isRC())
    {
        mod.seq = reverseComplement(mod.seq);
        mod.ovr.match.canonize();
    }

    // Initialize the offset value, the amount that a coordinate 
    // for the non-root sequence must be shifted so that
    // the sequences are aligned
    mod.offset = mod.ovr.match.inverseTranslate(0);
    m_overlaps.push_back(mod);
}

//
void MultiOverlap::add(const MOData& mod)
{
    assert(mod.ovr.id[0] == m_rootID);
    m_overlaps.push_back(mod);
}

//
void MultiOverlap::updateRootSeq(const std::string& newSeq)
{
    m_rootSeq = newSeq;
}

//
Overlap MultiOverlap::getOverlap(size_t idx) const
{
    assert(idx < m_overlaps.size());
    return m_overlaps[idx].ovr;
}



// Returns true if the overlap at idx has the include flag set
int MultiOverlap::getPartition(size_t idx) const
{
    assert(idx < m_overlaps.size());
    return m_overlaps[idx].partitionID;
}

// Returns true if the overlap at idx has the include flag set
void MultiOverlap::setPartition(size_t idx, int p)
{
    assert(idx < m_overlaps.size());
    m_overlaps[idx].partitionID = p;
}


//
bool MultiOverlap::qcCheck() const
{
    for(size_t i = 0; i < m_rootSeq.size(); ++i)
    {
        AlphaCount64 ac = getAlphaCount(i);
        size_t callSupport = ac.get(m_rootSeq[i]);
        if(callSupport < 2)
            return false;
    }
    return true;
}

//
void MultiOverlap::countOverlaps(size_t& prefix_count, size_t& suffix_count) const
{
    prefix_count = 0;
    suffix_count = 0;
    for(size_t i = 0; i < m_overlaps.size(); ++i)
    {
        if(!m_overlaps[i].ovr.match.coord[0].isContained())
        {
            if(m_overlaps[i].ovr.match.coord[0].isLeftExtreme())
                ++prefix_count;
            if(m_overlaps[i].ovr.match.coord[0].isRightExtreme())
                ++suffix_count;
        }
    }
}


//
int MultiOverlap::calculateCoverageOverlap()
{
    int rightmost_prefix = -1;
    int leftmost_suffix = m_rootSeq.length();

    for(size_t i = 0; i < m_overlaps.size(); ++i)
    {
        if(!m_overlaps[i].ovr.match.coord[0].isContained())
        {
            if(m_overlaps[i].ovr.match.coord[0].isLeftExtreme())
            {
                // Prefix overlap
                int end = m_overlaps[i].ovr.match.coord[0].interval.end;
                if(end > rightmost_prefix)
                    rightmost_prefix = end;
            }
            else
            {
                // Suffix overlap
                int start = m_overlaps[i].ovr.match.coord[0].interval.start;
                if(start < leftmost_suffix)
                    leftmost_suffix = start;
            }
        }
    }

    if(leftmost_suffix > rightmost_prefix)
        return 0;
    else
        return rightmost_prefix - leftmost_suffix + 1;
}

// Return the base in mod that matches the base at
// idx in the root seq. If mod does not overlap 
// the root at this position, returns '\0'
char MultiOverlap::getMODBase(const MOData& mod, int idx) const
{
    int trans_idx = idx - mod.offset;
    if(trans_idx >= 0 && size_t(trans_idx) < mod.seq.size())
        return mod.seq[trans_idx];
    else
        return '\0';
}

// Get an AlphaCount64 representing the nucleotides
// observed at the given column
AlphaCount64 MultiOverlap::getAlphaCount(int idx) const
{
    AlphaCount64 ac;

    ac.increment(m_rootSeq[idx]);
    for(size_t i = 0; i < m_overlaps.size(); ++i)
    {
        const MOData& curr = m_overlaps[i];
        // translate idx into the frame of the current sequence
        int trans_idx = idx - curr.offset;
        if(trans_idx >= 0 && size_t(trans_idx) < curr.seq.size())
        {
            ac.increment(curr.seq[trans_idx]);
        }
    }
    return ac;
}


// Print the MultiOverlap groups specified by the IntVec to stdout
void MultiOverlap::printGroups()
{
    for(int i = 0; i <= 1; ++i)
    {
        MultiOverlap mo_group(m_rootID, m_rootSeq);
        for(size_t j = 0; j < m_overlaps.size(); ++j)
        {
            const MOData& curr = m_overlaps[j];
            if(curr.partitionID == i)
                mo_group.add(curr);
        }
        std::cout << "MO GROUP " << i << "\n";
        mo_group.print();
    }
}

// Print the MultiOverlap to stdout
void MultiOverlap::print(int default_padding, int max_overhang)
{
    std::sort(m_overlaps.begin(), m_overlaps.end(), MOData::sortOffset);
    std::cout << "\nDrawing overlaps for read " << m_rootID << "\n";
    int root_len = int(m_rootSeq.size());
    
    // Print the root row at the bottom
    printRow(default_padding, max_overhang, root_len, 0, root_len, 0, 0.0f, m_rootSeq, m_rootID);

    for(size_t i = 0; i < m_overlaps.size(); ++i)
    {
        const MOData& curr = m_overlaps[i];
        int overlap_len = curr.ovr.match.getMaxOverlapLength();
        int nd = curr.ovr.match.countDifferences(m_rootSeq, curr.seq);
        double score = static_cast<double>(nd) / overlap_len;
        printRow(default_padding, max_overhang, root_len, 
                 curr.offset, overlap_len, nd, 
                 score, curr.seq, curr.ovr.id[1].c_str());
    }    
}

// Print a single row of a multi-overlap to stdout
void MultiOverlap::printRow(int default_padding, int max_overhang, 
                            int root_len, int offset, int overlap_len, int pid, 
                            double score, const std::string& seq, const std::string& id)
{
    int c_len = seq.length();

    // This string runs from c_offset to c_offset + len
    // Clip the string at -max_overhang to root_len + max_overhang
    int left_clip = std::max(offset, -max_overhang);
    int right_clip = std::min(offset + c_len, root_len + max_overhang);
    
    // translate the clipping coordinates to the string coords
    int t_left_clip = left_clip - offset;
    int t_right_clip = right_clip - offset;
    
    // Calculate the length of the left padding
    int padding = default_padding + left_clip;
    std::string leader = (t_left_clip > 0) ? "..." : "";
    std::string trailer = (t_right_clip < c_len) ? "..." : ""; 
    std::string clipped = seq.substr(t_left_clip, t_right_clip - t_left_clip);
    padding -= leader.size();

    assert(padding >= 0);
    std::string padding_str(padding, ' ');
    std::string outstr = padding_str + leader + clipped + trailer;
    printf("%s\t%d\t%d\t%lf\tID:%s\n", outstr.c_str(), overlap_len, pid, score, id.c_str());
}

void MultiOverlap::printMasked()
{
    printf("ROOT\t%s\t%s\n", m_rootSeq.c_str(), m_rootID.c_str());

    for(size_t j = 0; j < m_overlaps.size(); ++j)
    {
        std::string out;
        for(size_t i = 0; i < m_rootSeq.length(); ++i)
        {
            char b = getMODBase(m_overlaps[j], i);
            if(b == '\0')
                out.push_back('.');
            else if(b == m_rootSeq[i])
                out.push_back('=');
            else
                out.push_back(b);
        }
        printf("OVRL\t%s\t%s\n", out.c_str(), m_overlaps[j].ovr.id[1].c_str());
    }
}

//
bool MultiOverlap::MOData::sortOffset(const MOData& a, const MOData& b)
{
    return a.offset < b.offset;
}

// Sort by the ID of the non-root sequence, which is 
// guarenteed to bethe second id in the overlap structure
bool MultiOverlap::MOData::sortID(const MOData& a, const MOData& b)
{
    assert(a.ovr.id[0] == b.ovr.id[0]);
    return a.ovr.id[1] < b.ovr.id[1];
}

