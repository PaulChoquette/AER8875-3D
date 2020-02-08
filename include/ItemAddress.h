#pragma once

#include "CodeContract.h"
#include "StandardIntegralTypes.h"

namespace tecplot {

/**
 * Calculates the size of the identifier in bits.
 */
#define SIZE_IN_BITS(identifier) (sizeof(identifier)*8)

/**
 * Macro to verify that a item's partition has not exceeded the allowable number of bits
 * required to be represented by a tecplot::ItemAddress partition.
 */
#define VALID_ITEM_ADDRESS_PARTITION(itemAddressPartition) \
    ((itemAddressPartition) != tecplot::ItemAddress::INVALID_PARTITION && \
     (SIZE_IN_BITS(tecplot::ItemAddress::Partition_t) == tecplot::ItemAddress::PartitionBitSize || \
      (SIZE_IN_BITS(tecplot::ItemAddress::Partition_t) > tecplot::ItemAddress::PartitionBitSize && \
       uint64_t(itemAddressPartition) >> tecplot::ItemAddress::PartitionBitSize == uint64_t(0))))

/**
 * Macro to verify that a item's subzone offset has not exceeded the allowable number of bits
 * required to be represented by a tecplot::ItemAddress subzone offset.
 */
#define VALID_ITEM_ADDRESS_SUBZONE_OFFSET(itemAddressSubzoneOffset) \
    ((itemAddressSubzoneOffset) != tecplot::ItemAddress::INVALID_SUBZONE_OFFSET && \
     (SIZE_IN_BITS(tecplot::ItemAddress::SubzoneOffset_t) == tecplot::ItemAddress::SubzoneOffsetBitSize || \
      (SIZE_IN_BITS(tecplot::ItemAddress::SubzoneOffset_t) > tecplot::ItemAddress::SubzoneOffsetBitSize && \
       uint64_t(itemAddressSubzoneOffset) >> tecplot::ItemAddress::SubzoneOffsetBitSize == uint64_t(0))))

/**
 * Macro to verify that a subzone offset has not exceeded the allowable number of bits
 * required to be represented by a tecplot::ItemAddress item.
 */
#define VALID_ITEM_ADDRESS_ITEM_OFFSET(itemAddressItemOffset) \
    ((itemAddressItemOffset) != tecplot::ItemAddress::INVALID_ITEM_OFFSET && \
     (SIZE_IN_BITS(tecplot::ItemAddress::ItemOffset_t) == tecplot::ItemAddress::ItemOffsetBitSize || \
      (SIZE_IN_BITS(tecplot::ItemAddress::ItemOffset_t) > tecplot::ItemAddress::ItemOffsetBitSize && \
       uint64_t(itemAddressItemOffset) >> tecplot::ItemAddress::ItemOffsetBitSize == uint64_t(0))))

/**
 * Macro to verify that an item address is valid.
 */
#define VALID_UNIFORM_ITEM_ADDRESS(itemAddress) \
    ((itemAddress).addressType() == tecplot::ItemAddress::UniformAddressType)
#define VALID_SZL_ITEM_ADDRESS(itemAddress) \
    ((itemAddress).addressType() == tecplot::ItemAddress::SzlAddressType && \
     VALID_ITEM_ADDRESS_PARTITION((itemAddress).partition()) && \
     VALID_ITEM_ADDRESS_SUBZONE_OFFSET((itemAddress).subzoneOffset()) && \
     VALID_ITEM_ADDRESS_ITEM_OFFSET((itemAddress).itemOffset()))
#define VALID_ITEM_ADDRESS(itemAddress) \
    (VALID_UNIFORM_ITEM_ADDRESS(itemAddress) || \
     VALID_SZL_ITEM_ADDRESS(itemAddress))
#define VALID_ITEM_ADDRESS_SUBZONE_ADDRESS(subzoneAddress) \
    (VALID_ITEM_ADDRESS_PARTITION((subzoneAddress).partition()) && \
     VALID_ITEM_ADDRESS_SUBZONE_OFFSET((subzoneAddress).subzoneOffset()))


/**
 * An item address that support partition/subzoneOffset/item and uniform addressing. The class
 * represents the address of an item (cell, node or face) within a subzone using a partition and
 * subzone offset, or uniform addressing of a classic block of data. The copy performance of the
 * class is nearly identical to a word of the same size.
 */
class ItemAddress
{
public:
    /**
     * Type of the partition containing the subzone. All partition values must fit within the number
     * of bits specified by PartitionBitSize.
     */
    typedef uint32_t Partition_t;

    /**
     * For counting all subzones over all partitions, which could get over 4.29 billion subzones at
     * some point.
     */
    typedef uint64_t SubzoneIndex_t;

    /**
     * Type of the subzone containing the cell, node, or face item. All subzone offset values must
     * fit within the number of bits specified by SubzoneOffsetBitSize.
     */
    typedef uint32_t SubzoneOffset_t;

    /**
     * Type of the cell, node, or face item within the subzone. All item values must fit within the
     * number of bits specified by ItemOffsetBitSize.
     */
    typedef uint16_t ItemOffset_t;

    /**
     * The type used to hold a uniform address of an item. This allows non-SZL addressing to be
     * held by an ItemAddress for that existing algorithms in Tecplot work with both forms of data.
     * The uniform offset address is passed with a int64_t because it needs to hold uniform node
     * value addresses which are NumElements*PtsPerElement in size which is 2B*8.
     */
    typedef int64_t UniformOffset_t;

    /**
     * Number of bits needed to represent the integer identifying the addressing scheme used by
     * ItemAddress.
     */
    static uint32_t const AddressTypeBitSize = 1U;

    /**
     * Number of bits needed to represent the partition.
     */
    static uint32_t const PartitionBitSize = 23U;

    /**
     * Number of bits needed to represent the subzone offset.
     */
    static uint32_t const SubzoneOffsetBitSize = 32U;

    /**
     * Number of bits needed to represent the item within the subzone.
     */
    static uint32_t const ItemOffsetBitSize = 8U;

    /**
     * Number of bits needed to represent a uniform offset that doesn't use partition/subzone/item
     * (aka non-uniform) addressing. The uniform offset address is held with a 48 bits because it
     * needs to hold uniform node value addresses which are NumElements*PtsPerElement in size which
     * is 2B*8.
     */
    static uint32_t const UniformOffsetBitSize = 63U;

    static uint32_t const SzlAddressType     = 0U;
    static uint32_t const UniformAddressType = 1U;

    /**
     * Invalid partition value. Since we are using unsigned integers we sacrifice the last item.
     */
    static Partition_t const INVALID_PARTITION =
        static_cast<Partition_t>((uint64_t(1) << PartitionBitSize) - uint64_t(1));

    /**
     * To construct a valid ItemAddress when we neither know nor care what partition we're in.
     */
    static Partition_t const UNKNOWN_PARTITION = 0;

    /**
     * Maximum value for a item address partition. The maximum value is one less than the invalid index.
     */
    static Partition_t const MAX_PARTITION =
        INVALID_PARTITION - static_cast<Partition_t>(1);

    /**
     * Invalid subzone offset value. Since we are using unsigned integers we sacrifice the last item.
     */
    static SubzoneOffset_t const INVALID_SUBZONE_OFFSET =
        static_cast<SubzoneOffset_t>((uint64_t(1) << SubzoneOffsetBitSize) - uint64_t(1));

    /**
     * Maximum value for a item address subzone offset. The maximum value is one less than the invalid index.
     */
    static SubzoneOffset_t const MAX_SUBZONE_OFFSET =
        INVALID_SUBZONE_OFFSET - static_cast<SubzoneOffset_t>(1);

    /**
     * Invalid item offset value. Unlike the invalid partition or subzone offset which is one less
     * than the maximum value that can be represented by the bits, we don't have that flexibility
     * with item offsets because we need all 8 bits to represent 256 offsets so technically there is
     * no invalid 8 bit item offset however since the value is passed around as an ItemOffset_t
     * which is wider than 8 bits we can declare an invalid item, it just can be used to test an 8
     * bit item offset.
     */
    static ItemOffset_t const INVALID_ITEM_OFFSET =
        static_cast<ItemOffset_t>((uint64_t(1) << ItemOffsetBitSize));

    /**
     * Maximum value for a item address item. See INVALID_ITEM_OFFSET note.
     */
    static ItemOffset_t const MAX_ITEM_OFFSET =
        INVALID_ITEM_OFFSET - static_cast<ItemOffset_t>(1);

    /**
     */
    class SubzoneAddress
    {
    public:
        /**
         * Constructs an uninitialized subzone address. This is intentional so that a subzone
         * address has similar performance to a 64 bit word.
         * IMPORTANT:
         *     For performance reasons the member variables are not initialized when constructing using
         *     the default constructor, much in the way a 64 bit word is not initialized. This is
         *     particularly important when allocating large arrays of subzone addresses where the
         *     number of items being added to the array is not known upfront. In these cases the
         *     Linux/Windows memory allocator will not commit actual RAM until the memory is touched
         *     so we don't want to initialize items that may never be used.
         */
        inline SubzoneAddress()
        {
            INVARIANT(SIZE_IN_BITS(*this) == 64U);
        }

        inline SubzoneAddress(
            Partition_t     partition,
            SubzoneOffset_t subzoneOffset)
            : m_partition(partition)
            , m_subzoneOffset(subzoneOffset)
        {
            INVARIANT(SIZE_IN_BITS(*this) == 64U);
            REQUIRE(VALID_ITEM_ADDRESS_PARTITION(m_partition));
            REQUIRE(VALID_ITEM_ADDRESS_SUBZONE_OFFSET(m_subzoneOffset));
        }

        inline Partition_t partition() const
        {
            REQUIRE(VALID_ITEM_ADDRESS_PARTITION(m_partition));
            REQUIRE(VALID_ITEM_ADDRESS_SUBZONE_OFFSET(m_subzoneOffset));
            return m_partition;
        }

        inline SubzoneOffset_t subzoneOffset() const
        {
            REQUIRE(VALID_ITEM_ADDRESS_PARTITION(m_partition));
            REQUIRE(VALID_ITEM_ADDRESS_SUBZONE_OFFSET(m_subzoneOffset));
            return m_subzoneOffset;
        }

        inline bool operator==(SubzoneAddress const& other) const
        {
            REQUIRE(VALID_ITEM_ADDRESS_PARTITION(m_partition));
            REQUIRE(VALID_ITEM_ADDRESS_SUBZONE_OFFSET(m_subzoneOffset));
            return m_partition == other.m_partition && m_subzoneOffset == other.m_subzoneOffset;
        }

    private:
        Partition_t     m_partition;
        SubzoneOffset_t m_subzoneOffset;
    };

    /**
     * Constructs an uninitialized item address. This is intentional so that an item address has
     * similar performance to a 64 bit word.
     * IMPORTANT:
     *     For performance reasons the member variables are not initialized when constructing using
     *     the default constructor, much in the way a 64 bit word is not initialized. This is
     *     particularly important when allocating large arrays of item addresses where the number of
     *     items being added to the array is not known upfront. In these cases the Linux/Windows
     *     memory allocator will not commit actual RAM until the memory is touched so we don't want
     *     to initialize items that may never be used.
     */
    inline ItemAddress()
    {
        INVARIANT(validInvariants());
    }

    /**
     * Constructs an item address from a partition/subzone/item.
     * @param partition
     *     The partition containing the subzone. All partition values must fit within the number of
     *     bits specified by PartitionBitSize.
     * @param subzoneOffset
     *     The offset of the subzone containing the cell, node, or face item. All subzone offset
     *     values must fit within the number of bits specified by SubzoneOffsetBitSize.
     * @param itemOffset
     *     The cell, node, or face item within the subzone. All item values must fit within the
     *     number of bits specified by ItemOffsetBitSize.
     * @pre VALID_ITEM_ADDRESS_PARTITION(partition)
     * @pre VALID_ITEM_ADDRESS_SUBZONE_OFFSET(subzoneOffset)
     * @pre VALID_ITEM_ADDRESS_ITEM_OFFSET(itemOffset)
     */
    inline ItemAddress(
        Partition_t     partition,
        SubzoneOffset_t subzoneOffset,
        ItemOffset_t    itemOffset)
    {
        INVARIANT(validInvariants());
        REQUIRE(VALID_ITEM_ADDRESS_PARTITION(partition));
        REQUIRE(VALID_ITEM_ADDRESS_SUBZONE_OFFSET(subzoneOffset));
        REQUIRE(VALID_ITEM_ADDRESS_ITEM_OFFSET(itemOffset));

        /*
         * IMPORTANT:
         *     If the bit fields do not consume all 64 bits we would be forced to first
         *     initialize the entire 64 bit word to zero so that the equality test, see
         *     operator== will succeed: i.e.:
         *         _.rawBits = 0;
         */
        CHECK(AddressTypeBitSize + UniformOffsetBitSize == SIZE_IN_BITS(ItemAddress));
        CHECK(AddressTypeBitSize + PartitionBitSize + SubzoneOffsetBitSize + ItemOffsetBitSize == SIZE_IN_BITS(ItemAddress));

        m.addressType                  = SzlAddressType;
        m.szlItemAddress.partition     = partition;
        m.szlItemAddress.subzoneOffset = subzoneOffset;
        m.szlItemAddress.itemOffset    = itemOffset;
    }

    /**
     * Constructs an item address from a subzone address and an node, cell, or face within the subzone.
     * @param subzoneAddress
     *     subzone address in which the item is located
     * @param itemOffset
     *     item offset of the node, cell or face within the subzone
     */
    inline ItemAddress(
        SubzoneAddress subzoneAddress,
        ItemOffset_t   itemOffset)
    {
        INVARIANT(validInvariants());
        REQUIRE(VALID_ITEM_ADDRESS_SUBZONE_ADDRESS(subzoneAddress));
        REQUIRE(VALID_ITEM_ADDRESS_ITEM_OFFSET(itemOffset));

        /*
         * IMPORTANT:
         *     If the bit fields do not consume all 64 bits we would be forced to first
         *     initialize the entire 64 bit word to zero so that the equality test, see
         *     operator== will succeed: i.e.:
         *         m.rawBits = 0;
         */
        CHECK(AddressTypeBitSize + UniformOffsetBitSize == SIZE_IN_BITS(ItemAddress));
        CHECK(AddressTypeBitSize + PartitionBitSize + SubzoneOffsetBitSize + ItemOffsetBitSize == SIZE_IN_BITS(ItemAddress));

        m.addressType                  = SzlAddressType;
        m.szlItemAddress.partition     = subzoneAddress.partition();
        m.szlItemAddress.subzoneOffset = subzoneAddress.subzoneOffset();
        m.szlItemAddress.itemOffset    = itemOffset;
    }

    /**
     * Constructs a item address from a uniform offset.
     */
    inline explicit ItemAddress(UniformOffset_t uniformOffset)
    {
        INVARIANT(validInvariants());

        /*
         * IMPORTANT:
         *     If the bit fields do not consume all 64 bits we would be forced to first
         *     initialize the entire 64 bit word to zero so that the equality test, see
         *     operator== will succeed: i.e.:
         *         m.rawBits = 0;
         */
        CHECK(AddressTypeBitSize + UniformOffsetBitSize == SIZE_IN_BITS(ItemAddress));
        CHECK(AddressTypeBitSize + PartitionBitSize + SubzoneOffsetBitSize + ItemOffsetBitSize == SIZE_IN_BITS(ItemAddress));

        m.addressType           = UniformAddressType;
        m.uniformAddress.offset = uniformOffset;
    }

    /**
     * @return
     *     true if the item address is using uniform addressing, false otherwise
     */
    inline bool isUniform() const
    {
        ENSURE(m.addressType == SzlAddressType || m.addressType == UniformAddressType);
        return m.addressType == UniformAddressType;
    }

    /**
     * @return
     *     true if the item address is using single-item addressing, false otherwise
     */
    inline bool isSzlItem() const
    {
        ENSURE(m.addressType == SzlAddressType || m.addressType == UniformAddressType);
        return m.addressType == SzlAddressType;
    }

    /**
     * @return
     *     The type of addressing of the item address
     */
    inline uint32_t addressType() const
    {
        ENSURE(m.addressType == SzlAddressType || m.addressType == UniformAddressType);
        return m.addressType;
    }

    /**
     * @return
     *     The subzone address of this item address.
     */
    inline SubzoneAddress subzoneAddress() const
    {
        REQUIRE(m.addressType == SzlAddressType);
        return SubzoneAddress(m.szlItemAddress.partition, m.szlItemAddress.subzoneOffset);
    }

    /**
     * @return
     *     partition in which the item address is located
     */
    inline Partition_t partition() const
    {
        REQUIRE(m.addressType == SzlAddressType);
        ENSURE(VALID_ITEM_ADDRESS_PARTITION(m.szlItemAddress.partition));
        return m.szlItemAddress.partition;
    }

    /**
     * @return
     *     subzone in which the item address is located
     */
    inline SubzoneOffset_t subzoneOffset() const
    {
        REQUIRE(m.addressType == SzlAddressType);
        ENSURE(VALID_ITEM_ADDRESS_SUBZONE_OFFSET(m.szlItemAddress.subzoneOffset));
        return m.szlItemAddress.subzoneOffset;
    }

    /**
     * @return
     *     node, cell, or face within the subzone
     */
    inline ItemOffset_t itemOffset() const
    {
        REQUIRE(m.addressType == SzlAddressType);
        ENSURE(VALID_ITEM_ADDRESS_ITEM_OFFSET(m.szlItemAddress.itemOffset));
        return m.szlItemAddress.itemOffset;
    }

    /**
     * @return
     *     uniform offset of the item (cell, node, or face)
     */
    inline UniformOffset_t uniformOffset() const
    {
        REQUIRE(isUniform());
        return m.uniformAddress.offset;
    }

    /**
     * @return
     *     Item address raw bit pattern. Note that this is not the same as a uniform offset. In
     *     addition it is important to note that the raw bit pattern varies between compilers and
     *     should not be stored in a file for later use. The raw bits are exposed primarily so that
     *     an item addresses can be temporarily stored in a 64 bit field
     */
    inline uint64_t toRawBits() const
    {
        return m.rawBits;
    }

    /**
     * @param
     *     Raw bits of an item address.
     * @return
     *     An item address constructed from the raw bits of another item address. Note that this is
     *     not the same as a uniform offset. In addition it is important to note that the raw bit
     *     pattern varies between compilers and should not be stored in a file for later use. The
     *     raw bits are exposed primarily so that an item addresses can be temporarily stored in a
     *     64 bit field.
     */
    static ItemAddress fromRawBits(uint64_t rawBits)
    {
        ItemAddress result;
        result.m.rawBits = rawBits;
        ENSURE(VALID_ITEM_ADDRESS(result));
        return result;
    }

    /**
     * See binary comparison operations for other operations.
     * @return
     *     true if the two items have bitwise equality
     */
    inline bool operator==(ItemAddress const& other) const
    {
        REQUIRE(VALID_ITEM_ADDRESS(other));
        return m.rawBits == other.m.rawBits;
    }

private:
    #if !defined NO_ASSERTS
    static bool validInvariants()
    {
        return (SIZE_IN_BITS(ItemAddress) <= 64U && // ...otherwise consider passing ItemAddress' by reference instead of by value
                SIZE_IN_BITS(Partition_t)     >= PartitionBitSize     &&
                SIZE_IN_BITS(SubzoneOffset_t) >= SubzoneOffsetBitSize &&
                SIZE_IN_BITS(ItemOffset_t)    >  ItemOffsetBitSize    && // ...unlike the other bit fields we need all 8 bits so ItemOffset_t must be larger for correct loop termination
                SIZE_IN_BITS(UniformOffset_t) >= UniformOffsetBitSize &&
                AddressTypeBitSize + UniformOffsetBitSize == SIZE_IN_BITS(ItemAddress) &&
                AddressTypeBitSize + PartitionBitSize + SubzoneOffsetBitSize + ItemOffsetBitSize == SIZE_IN_BITS(ItemAddress));
    }
    #endif

    /**
     * IMPORTANT:
     *     Courtesy of Visual Studio, all bit field members within a structure must be of the same
     *     type to instruct the compiler to pack into one 64 bit space. Because of this limitation,
     *     and because we must preserve the sign bit required by uniform addresses in Tecplot (e.g.
     *     we use -1 and -2 as flags in cell and node numbering), we have to make them all members
     *     of the UniformAddress_s structure signed 64 bit integers. The same is true for the
     *     SzlItemAddress_s structure except all the members must be 64 bit unsigned values.
     *     Basically, a structure containing bit fields must either all be signed or unsigned to
     *     work. If this were not the case we would have to abandon bit fields altogether and
     *     perform the bit twiddling ourselves. After all the grief of supposedly portable bit
     *     fields I'm tempted to abandon them anyway.
     */

    struct SzlItemAddress_s
    {
        // IMPORTANT: see note above for explanation of type choices for bit fields.
        uint64_t addressType:AddressTypeBitSize;
        uint64_t partition:PartitionBitSize;
        uint64_t subzoneOffset:SubzoneOffsetBitSize;
        uint64_t itemOffset:ItemOffsetBitSize;
    };

    struct UniformAddress_s
    {
        // IMPORTANT: see note above for explanation of type choices for bit fields.
        int64_t addressType:AddressTypeBitSize;
        int64_t offset:UniformOffsetBitSize;
    };

    union
    {
        // IMPORTANT: see note above for explanation of type choices for bit fields.
        uint64_t         addressType:AddressTypeBitSize;
        SzlItemAddress_s szlItemAddress;
        UniformAddress_s uniformAddress;
        uint64_t         rawBits; // ...used for bitwise equality
    } m;
};


/*
 * Binary infix comparison operators for ItemAddress. ItemAddress supplies a bitwise unary
 * operator== and a binary operator<. All other comparison operators are written in terms of these
 * two operators.
 */

inline bool operator<(
    ItemAddress const& lhs,
    ItemAddress const& rhs)
{
    REQUIRE(VALID_ITEM_ADDRESS(lhs));
    REQUIRE(VALID_ITEM_ADDRESS(rhs));

    uint32_t const lhsType = lhs.addressType();
    uint32_t const rhsType = rhs.addressType();
    if (lhsType != rhsType)
    {
        return lhsType < rhsType;
    }
    else if (lhsType == tecplot::ItemAddress::UniformAddressType)
    {
        return lhs.uniformOffset() < rhs.uniformOffset();
    }
    else if (lhs.partition() == rhs.partition())
    {
        if (lhs.subzoneOffset() == rhs.subzoneOffset())
            return lhs.itemOffset() < rhs.itemOffset();
        else
            return lhs.subzoneOffset() < rhs.subzoneOffset();
    }
    else
    {
        return lhs.partition() < rhs.partition();
    }
}

inline bool operator!=(
    ItemAddress const& lhs,
    ItemAddress const& rhs)
{
    REQUIRE(VALID_ITEM_ADDRESS(lhs));
    REQUIRE(VALID_ITEM_ADDRESS(rhs));
    return !(lhs == rhs);
}

inline bool operator<=(
    ItemAddress const& lhs,
    ItemAddress const& rhs)
{
    REQUIRE(VALID_ITEM_ADDRESS(lhs));
    REQUIRE(VALID_ITEM_ADDRESS(rhs));
    return lhs < rhs || lhs == rhs;
}

inline bool operator>(
    ItemAddress const& lhs,
    ItemAddress const& rhs)
{
    REQUIRE(VALID_ITEM_ADDRESS(lhs));
    REQUIRE(VALID_ITEM_ADDRESS(rhs));
    return !(lhs < rhs) && !(lhs == rhs);
}

inline bool operator>=(
    ItemAddress const& lhs,
    ItemAddress const& rhs)
{
    REQUIRE(VALID_ITEM_ADDRESS(lhs));
    REQUIRE(VALID_ITEM_ADDRESS(rhs));
    return !(lhs < rhs);
}

inline bool operator<(
    ItemAddress::SubzoneAddress const& lhs,
    ItemAddress::SubzoneAddress const& rhs)
{
    REQUIRE(VALID_ITEM_ADDRESS_SUBZONE_ADDRESS(lhs));
    REQUIRE(VALID_ITEM_ADDRESS_SUBZONE_ADDRESS(rhs));

    if (lhs.partition() == rhs.partition())
        return lhs.subzoneOffset() < rhs.subzoneOffset();
    else
        return lhs.partition() < rhs.partition();
}

inline bool operator!=(
    ItemAddress::SubzoneAddress const& lhs,
    ItemAddress::SubzoneAddress const& rhs)
{
    REQUIRE(VALID_ITEM_ADDRESS_SUBZONE_ADDRESS(lhs));
    REQUIRE(VALID_ITEM_ADDRESS_SUBZONE_ADDRESS(rhs));
    return !(lhs == rhs);
}

inline bool operator<=(
    ItemAddress::SubzoneAddress const& lhs,
    ItemAddress::SubzoneAddress const& rhs)
{
    REQUIRE(VALID_ITEM_ADDRESS_SUBZONE_ADDRESS(lhs));
    REQUIRE(VALID_ITEM_ADDRESS_SUBZONE_ADDRESS(rhs));
    return lhs < rhs || lhs == rhs;
}

inline bool operator>(
    ItemAddress::SubzoneAddress const& lhs,
    ItemAddress::SubzoneAddress const& rhs)
{
    REQUIRE(VALID_ITEM_ADDRESS_SUBZONE_ADDRESS(lhs));
    REQUIRE(VALID_ITEM_ADDRESS_SUBZONE_ADDRESS(rhs));
    return !(lhs < rhs) && !(lhs == rhs);
}

inline bool operator>=(
    ItemAddress::SubzoneAddress const& lhs,
    ItemAddress::SubzoneAddress const& rhs)
{
    REQUIRE(VALID_ITEM_ADDRESS_SUBZONE_ADDRESS(lhs));
    REQUIRE(VALID_ITEM_ADDRESS_SUBZONE_ADDRESS(rhs));
    return !(lhs < rhs);
}

}
