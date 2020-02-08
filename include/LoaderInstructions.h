#pragma once

#include "ThirdPartyHeadersBegin.h"
    #include <memory>
    #include <string>
    #include <vector>
#include "ThirdPartyHeadersEnd.h"

#include "ClassMacros.h"
#include "StringList.h"
#include "TecUtil.h"
#include "tecutildataio_utils_Exports.h"


namespace tecplot { namespace tecutildataio_utils {

/**
 * Holds onto the instructions while setting up the callbacks.
 */
class tecutildataio_utils_API LoaderInstructions
{
public:
    /**
     * @throws std::bad_alloc if there are insufficient resources.
     * @throws InvalidInstructionsException if instructions do not match the specification
     */
    LoaderInstructions(
        std::string const& dataLoaderName,
        StringList_pa      instructions);

    /**
     * @throws std::bad_alloc if there are insufficient resources.
     */
    explicit LoaderInstructions(std::vector<std::shared_ptr<std::string const>> const& filenames);

    /**
     * @throws std::bad_alloc if there are insufficient resources
     */
    tecplot::toolbox::StringList getStringList() const;

    /**
     * @throws std::bad_alloc if there are insufficient resources.
     */
    std::vector<std::shared_ptr<std::string const>> const& filenames() const;

private:
    UNCOPYABLE_CLASS(LoaderInstructions);
    std::vector<std::shared_ptr<std::string const>> m_filenames;
};

}}
