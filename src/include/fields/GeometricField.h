#ifndef GEOMETRICFIELD_H
#define GEOMETRICFIELD_H

#include <vector>

#include "utilities/Reference.h"
#include "utilities/Utilities.h"
#include "iomanagment/InfoStream.h"
#include "iomanagment/RegistryFile.h"
#include "iomanagment/RegistryObject.h"
#include "mesh/Mesh.h"
#include "time/TimeChangeListner.h"
#include "time/Time.h"
#include "fields/PatchField.h"
#include "fields/PatchSelection.h"
#include "components/CmpTraits.h"
#include "utilities/Utilities.h"
#include "mesh/Mesh.h"
#include "ElementFieldBase.h"
#include "ContinousField.h"

namespace SEM { namespace field {
    
    /** ******************************************
     * \class GeometricField
     * Generic class for any entity filed, assocciated with mesh.
     * Supports automatic file read/write and control
     * over time change-store old results.
     * It's generic class that can be used in pde, which
     * keeps information abut boundary conditions applied
     * in file.
     ********************************************/
    //-----------------------------------------------------------------
    //  TODO - move timeChangeListner implementation to some higher level
    //         class, because initialize some member with NULL is not best
    //         solution
    //  TODO -perhaps this class shall have 2nd template argument to
    //        provide compile time chosing number of time to be stored.
    //  TODO - add "clone" method to Patch field and use it inside operator |=
    //-----------------------------------------------------------------
    template <typename T>
    class GeometricField :  public ContinousField<T>, public time::TimeChangeListner
    {
        REFERENCE_TYPE(GeometricField<T> )
    public:
        typedef PatchField<T> Patch;
        typedef std::vector<Patch*> BoundaryField;
        
    private:
        typedef numArray<T> Storage;
        typedef ContinousField<T> FieldBase;
        
        static int NCached;
        
        ///  \brief m_oldTimes - ordered old field acording to saving
        ///   (old values are stored as pointers. When time is changed,
        ///    then last pointer became the first one and rest pointers
        ///    are moved back. In this way only 1 assigment between field
        ///    values is made - from this class object to last(after reordering
        ///    to first) objcet pointed by pointer)
        
        
        /// \var m_oldTimes - orderd old field values.
        /// Old values are stored as pointers, to allow fast reordering. 
        /// When cacheCurrentField is called, pointers in std::vector 
        /// moves it's location one back. At first index=0 is stroed
        /// copy of field values for actual time, because field values
        /// may change at the same time iteration, and current time 
        /// values are assumed to be that which coresponds to field
        /// values at the moment when cacheCurrentField is called.
        /// Eg. cacheField(0)-cached values for actual time,
        /// cacheFiled(1) - cached values for previous time, and so on...
        std::vector<Storage*> m_oldTimes;
        
        
        /// \var m_boundaryFields:  collection of derived PatchFields
        BoundaryField m_boundaryFields;
        
        /// \var m_timeRegister:  object which controls caching current values,
        ///  need to be stored to allow unregistering when destructor would be
        ///  called.
        Time* m_timeRegister;
        
        /// \var m_mesh: reference to mesh object
        const mesh::Mesh & m_mesh;
        
    public:
        using iomanagment::RegistryObject::setRegistryFile;
        using FieldBase::slice;
        
        //------------------------------------------------------------------------------//
        //                      CONSTRUCTORS
        //------------------------------------------------------------------------------//
        /** ***************************************************************************
         * \brief GeometricField - create field basing on geometrical mesh size. Field
         *                         is not associated in any register which would control
         *                         read/write and field values caching operation
         * \param mesh           - geometrical object
         *******************************************************************************/
        GeometricField(const mesh::Mesh& mesh, const std::string& name = "GeometricField")
        : FieldBase(mesh.nodesNumber(),name), m_timeRegister(NULL), m_mesh(mesh)
        {
        }
        
        /** ******************************************************************************
         * \function GeometricField - construct field which is registred as time change listner
         *   (store filed result in cached filed for each time when "time" object would 
         *   decide-normaly each time iteration).
         *   Read/write operation are not automaticly performed.
         *  \param time - object which will trigger storing 
         *******************************************************************************/
        GeometricField(Time& time, const mesh::Mesh& mesh, const std::string& name = "GeometricField")
        :FieldBase(mesh.nodesNumber(),name), m_timeRegister(&time), m_mesh(mesh)
        {
            time.registerTimeListner(this);
        }
        
        /** ****************************************************************************** 
         * \brief GeometricField - construct geometirc field, which is registered only in RegistryFile
         *                         (object controling read/write with to specified file-files)
         * \param regFile        - associated file from/to wich data shall be read/write, according
         *                         to regFile specification. regFile is shared pointer. Use this
         *                         constructor if one regFile need to be shared between multiple
         *                         objects.
         * \param mesh           - geometrical object
         *********************************************************************************/
        explicit GeometricField(iomanagment::RegistryFile::ref regFile, const mesh::Mesh &mesh)
        : FieldBase(regFile, mesh.nodesNumber()), m_timeRegister(NULL), m_mesh(mesh)
        {
        }
        
        /** ****************************************************************************** 
         * \brief GeometricField - construct geometirc field, which is registered in RegistryFile
         *                         (object controling read/write with to specified file-files)
         *                         and as time change listner.(store field result in cached field
         *                         for each time when "time" object would decide-normaly each time
         *                         iteration)
         * \param time           - object which will trigger storing this field values as cached-field.
         *                         normaly it's time wich triggers it each timestep change
         * \param regFile        - associated file from/to wich data shall be read/write, according
         *                         to regFile specification. regFile is shared pointer. Use this
         *                         constructor if one regFile need to be shared between multiple
         *                         objects.
         * \param mesh           - geometrical object
         *********************************************************************************/
        explicit GeometricField(Time& time,iomanagment::RegistryFile::ref regFile, const mesh::Mesh &mesh)
        : FieldBase(regFile,mesh.nodesNumber()), m_timeRegister(&time), m_mesh(mesh)
        {
            time.registerTimeListner(this);
        }
        
        /** *******************************************************************************
         * \brief GeometricField - construct geometirc field, which is registered only in RegistryFile
         *                         (object controling read/write with to specified file-files)
         * \param time           - object which will trigger storing this field values as cached-field.
         *                         normaly it's time wich triggers it each timestep change
         * \param regFile        - associated file from/to wich data shall be read/write, according
         *                         to regFile specification.
         * \param mesh           - geometrical object
         * --------------------------------------------------------------------------------------
         * WARNING!!! this constructor for safety purpose shall not exist.
         * It shall be only used when regFilePtr will be passed to this instance.
         * It's intorduced for case :
         *         GeometricField<T>(new RegistryFile(...),mesh)--->much more readable
         * This constructor wraps raw pointer with boost::shared_ptr, so when
         * wrapping obj. will be destroyied then regFilePtr will be delete.
         *********************************************************************************/
        GeometricField(iomanagment::RegistryFile* regFilePtr, const mesh::Mesh &mesh)
        : FieldBase(mesh.nodesNumber(),regFilePtr->fileName()), m_timeRegister(NULL), m_mesh(mesh)
        {
            iomanagment::RegistryFile::ref regFile(regFilePtr);
            setRegistryFile(regFile);
        }
        
        /** *******************************************************************************
         * \brief GeometricField - construct geometirc field, which is registered in RegistryFile
         *                         (object controling read/write with to specified file-files)
         *                         and as time change listner.(store field result in cached field
         *                         for each time when "time" object would decide-normaly each time
         *                         iteration)
         * \param time           - object which will trigger storing this field values as cached-field.
         *                         normaly it's time wich triggers it each timestep change
         * \param regFile        - associated file from/to wich data shall be read/write, according
         *                         to regFile specification.
         * \param mesh           - geometrical object
         * --------------------------------------------------------------------------------------
         * WARNING!!! this constructor for safety purpose shall not exist.
         * It shall be only used when regFilePtr will be passed to this instance.
         * It's intorduced for case :
         *         GeometricField<T>(new RegistryFile(...),mesh)--->much more readable
         * This constructor wraps raw pointer with boost::shared_ptr, so when
         * wrapping obj. will be destroyied then regFilePtr will be delete.
         * *********************************************************************************/
        GeometricField(Time &time, iomanagment::RegistryFile* regFilePtr, const mesh::Mesh &mesh)
        : FieldBase(mesh.nodesNumber(),regFilePtr->fileName()), m_timeRegister(&time), m_mesh(mesh)
        {
            iomanagment::RegistryFile::ref regFile(regFilePtr);
            setRegistryFile(regFile);
            time.registerTimeListner(this);
        }
        
        /** *********************************************************************************
         * \brief GeometricField "easy" constructor
         *  Takes time variable from static Case object as object registry. The same variable
         *  is also used for time register. As default new FileRegistry is created, so
         *  user need to pay attension to not duplicate FileRegistry pointing to the same file on
         *  disc space(it wont be error but unnecessery overhead). As local path file would be
         *  used time-time folders on disc space.
         * \param fieldName- name for this field-the same as file on disc space
         * \param m - mesh object
         * \param read - file reading option
         * \param write - file writing option
         *********************************************************************************/
        GeometricField(const std::string &fieldName, const mesh::Mesh &m, iomanagment::READ read=iomanagment::READ_ONCE, iomanagment::WRITE write=iomanagment::AUTO)
        : FieldBase(m.nodesNumber(),fieldName),m_timeRegister(&Case::time()), m_mesh(m)
        {
            iomanagment::RegistryFile* regFilePtr = new iomanagment::RegistryFile
            (
                Case::time(),
                fieldName,
                Case::time().localPath(),
                read,
                write
            );
            iomanagment::RegistryFile::ref regFile(regFilePtr);
            setRegistryFile(regFile);
            Case::time().registerTimeListner(this);
        }
        
        /** *********************************************************************************
         * \brief GeometricField copy-constructor
         * \param other              - other the same type object
         *********************************************************************************/
        GeometricField(const GeometricField<T>& other)
        :FieldBase(other.mesh().nodesNumber(),other.name()), m_timeRegister(other.m_timeRegister), m_mesh(other.mesh())
        {
            setRegistryFile(other.registryFile());
            m_timeRegister->registerTimeListner(this);
        }
        
        //------------------------------------------------------------------------------//
        //                      DESTRUCTOR                                  //
        //------------------------------------------------------------------------------//
        /** *******************************************************************************
         * \brief ~GeometricField - destructor
         *                         clear old times and patches allocated on stack,
         *                         and unregister from time(unregistred from RegFile in
         *                         regObject destructor).
         *********************************************************************************/
        virtual ~GeometricField()
        {
            //delete Patches on stack
            deleteAll(m_boundaryFields);
            m_boundaryFields.clear();
            
            // delete cached values on stack
            deleteAll(m_oldTimes);
            m_oldTimes.clear();
            
            // unregister from time if this object was constructed as listner.
            if(m_timeRegister)
            {
                m_timeRegister->unregisterTimeListner(this);
            }
        }
        
        //------------------------------------------------------------------------------//
        //                      ASSIGMENT                                     //
        //------------------------------------------------------------------------------//
        /** ******************************************************************************
         * \brief operator = assign only internal field
         * \param other - other geometric field
         * \return     - this object reference
         * *********************************************************************************/
        GeometricField<T> & operator =(const GeometricField<T>& other)
        {
            FieldBase::operator=(other);
            return *this;
        }
        
        GeometricField<T> & operator =(const Storage& other)
        {
            FieldBase::operator=(other);
            return *this;
        }
        
//         /** *******************************************************************************
//          * \brief operator |= assign internal field and patches as well
//          * \param other - other geometirc field
//          * \return reference to this object
//          *********************************************************************************/
//         GeometricField<T> & operator |=(const GeometricField<T>& other)
//         {
//             FieldBase::operator=(other);
//             m_boundaryFields = other.m_boundaryFields;
//             return *this;
//         }
        
        /** *******************************************************************************
         * \brief operator |= assign selected patch type to all boundaries
         * \param patchType - name of patch type
         * \return reference to this object
         *********************************************************************************/
        GeometricField<T> & operator |= (const std::string & patchType)
        {
            deleteAll(m_boundaryFields);
            m_boundaryFields.clear();
            
            PatchFactory<T> patches;
            
            for(const mesh::Boundary & boundary : mesh().boundaryMesh())
            {
                m_boundaryFields.push_back( patches.select(patchType,*m_timeRegister,boundary) );
            }
        }
        
        
        /** *******************************************************************************
         * \brief operator |= assign value patchesValue to all patches 
         * \param patchesValue - value to assign
         * \return reference to this object
         *********************************************************************************/
        GeometricField<T> & operator |= (const T& patchesValue)
        {
            typename std::vector<Patch *>::iterator itr=m_boundaryFields.begin();
            for(; itr!=m_boundaryFields.end(); ++itr )   
            {
                *(*itr) = patchesValue;
            }
        }
        
        //------------------------------------------------------------------------------//
        //                      UTILITIES                                              
        //------------------------------------------------------------------------------//
        /** ******************************************************************************
         * \brief mesh - const mesh object getter
         * \return     - mesh object associated with this filed
         * *******************************************************************************/
        const mesh::Mesh & mesh() const { return m_mesh;}
        
        /** ******************************************************************************
         * \brief element - local element field values getter
         * \return     - numArray view mapped by node indexes
         * *******************************************************************************/
        numArrayIndexMapped<numArray<T> > element(size_t e) 
        {
            return slice(mesh()[e].indexVectorMask());
        }
        
        /** ******************************************************************************
         * \brief element - const local element field values getter
         * \return     - numArray view mapped by node indexes
         * *******************************************************************************/
        numArrayIndexMapped<const numArray<T> > element(size_t e) const
        {
            return slice(mesh()[e].indexVectorMask());
        }
        
        /** ******************************************************************************
         * boundrayFields getter
         * \return collection of patches
         * *******************************************************************************/
        BoundaryField & boundaryFields(){ return m_boundaryFields; }
        
        /** ******************************************************************************
         * boundrayFields const getter
         * \return const collection of patches
         * *******************************************************************************/
        const BoundaryField & boundaryFields() const { return m_boundaryFields;}
        
        /** ******************************************************************************
         * \brief containPatch - check if patch with specified name exist here
         * \param patchName    - name of boundary to check
         * \return             true/false if patch exist
         * *******************************************************************************/
        const bool containPatch(const std::string &patchName) const
        {
            for(Patch* p: m_boundaryFields)
            {
                if(p->boundaryEdges().name()==patchName)
                    return true;
            }
            return false;
        }
        
        /** ******************************************************************************
         * \brief patch const patch getter
         * \param patchName - boundary name for patch to get
         * \return         - patch field
         * *******************************************************************************/
        const Patch & patch(const std::string &patchName) const
        {
            using namespace iomanagment;
            
            for (Patch* p : m_boundaryFields)
            {
                if(p->boundaryEdges().name()==patchName)
                    return *p;
            }
            ErrorInFunction<<"There is no patch for boundary named "<<patchName<<endProgram;
        }
        
        Patch & patch(const std::string &patchName)
        {
            using namespace iomanagment;
            
            for (Patch* p : m_boundaryFields)
            {
                if(p->boundaryEdges().name()==patchName)
                    return *p;
            }
            ErrorInFunction<<"There is no patch for boundary named "<<patchName<<endProgram;
        }
        
        
        /** ******************************************************************************
         * \brief isCachedAll - checks if number of fields which shall be cached is already
         *                    achived.
         * *******************************************************************************/
        bool isCachedAll() const { return m_oldTimes.size() == NCached; }
        
        unsigned int oldFieldsNumber() const { return m_oldTimes.size(); }
        
        
        /** ******************************************************************************
         * \brief cacheCurrentField - triggers storing values of this filed in cache
         *                            last field in cache memory will be replaced by ahead
         *                            and first one in cache will be assigned by this field values
         *                            (reasigning is done on 1 field, in other fields pointers only
         *                             change reference index location).
         * *******************************************************************************/
        void cacheCurrentField()
        {
            if(isCachedAll())
            {
                //realocate last element to first and then assign to its values from "this" field
                Storage* tmp = m_oldTimes.back();      //take last
                m_oldTimes.pop_back();                    //remove last
                (*tmp) = *this;                           //reassign values to last
                m_oldTimes.insert(m_oldTimes.begin(),tmp);//insert last as first
            }
            else
            {
                // Set current cached field as first in old times
                m_oldTimes.insert(m_oldTimes.begin(), new Storage(*this));
            }
        }
        
        /** ******************************************************************************
         * \brief cachedField - const cached field values
         * \param prevIndex   - specify which prev values you want to access.
         *                      (now implementation allows indexes 0,1 values (prev and prev preve).
         * \return            - cached field.
         * *******************************************************************************/
        const Storage& oldField(int prevIndex) const
        {
            using namespace iomanagment;
            
            if(prevIndex>NCached-1)
                ErrorInFunction<<"There is no more then "<<NCached<<" chacedFields in GeometricField\n"
                <<"Can't return previous cached field at index "<<prevIndex
                <<endProgram;
            
            //If yet not stored enough fileds then return the oldest one
            if(prevIndex+1>m_oldTimes.size())
                prevIndex=m_oldTimes.size()-1;
            
            //if none stored yet return this field itself
            if( prevIndex < 0)
                return *this;
            else
                return *m_oldTimes[prevIndex];
        }
        
        /** ******************************************************************************
         * \brief relax - function to relax current values in field by some factor.
         *                values in this field ar set as result of this=prev+(this-prev)*facor.
         * \param factor- relaxation factor
         * *******************************************************************************/
        void relax(Scalar factor)
        {
            (*this) = oldField(0) + ( (*this) - oldField(0))*factor;
        }
        
        /** ******************************************************************************
         * \brief updateBoundaryConditions
         * function which calls update metod in all patches. Shall be called if 
         * boundary conditions values depends on some other field, which values has changed.
         * \note in future it shall be moved to be some kind of listner triggerd whenever
         * operator= was used with field object. No it must be manualy triggered
         * *******************************************************************************/
        void updateBoundaryConditions()
        {
            for(auto & patch: m_boundaryFields)
            {
                patch->updateValues();
            }
        }
        
        
        //----------------------------------------------------------------------//
        //          IO-MANAGMENT
        // NOTE: in future, patch selection shall be made by some sort of
        //       runtime selection map - holding key-instanceProvider
        //----------------------------------------------------------------------//
        /** ******************************************************************************
         * \brief read - RegistryObject implementation. Field is setuped by provided
         *               dictionary with data(aimed for automatic reading: file->dict->object)
         * \param dict - provided dictionaray with data
         * *******************************************************************************/
        virtual void read(const iomanagment::Dictionary &dict)
        {
            using namespace iomanagment;
            
            //read interanl field values
            dict.entry("internalField")>>*this;
            
            if(dict.hasSubDictionary("boundary"))
            {
                Dictionary& boundaries = dict.subDictionary("boundary");
                PatchFactory<T> patches;
                
                for(auto mapEntry : boundaries.subDictionaries())
                {
                    const Dictionary & boundary = *(mapEntry.second);
                    const DictEntry& type = boundary.entry("type");
                    
                    bool create=true;
                    
                    // if patch already exist and is diffrent type then remove old one
                    // if patch already exist and type is equal then update from file.
                    if(containPatch(boundary.name()))
                    {
                        typename BoundaryField::iterator it=m_boundaryFields.begin();
                        for(;it!=m_boundaryFields.end(); ++it)
                        {
                            if( (*it)->boundaryEdges().name()==boundary.name() )
                            {
                                if(type != (*it)->typeName()) //remove patch to assign new type
                                {
                                    m_boundaryFields.erase(it);
                                    delete *it;
                                }
                                else // only update patch
                                {
                                    create = false;
                                    (*it)->read(dict);
                                }
                                break;
                            }
                        }
                    }
                    
                    if(create)
                    {
                        m_boundaryFields.push_back(patches.select(type.value(),*m_timeRegister,mesh().boundaryAt(boundary.name() ) ) );
                        m_boundaryFields.back()->read(dict);
                    }
                }
            }
        }
        
        
        //note: patches are not registred so they will not also automaticly write theier data
        //      (in fact they can be registered with the same file as this geometric filed,
        //      becasuse there is the place where they are defined. But if they would be registred, then
        //      there will be a little bit more complicated work to print for console this file dictionary,
        //      while in  this way it will require only putting some dict into write method. It would be
        //      more complicated because in below method, there will be no calling write on patches, but
        //      this call would be managed by regFile).
        /** ******************************************************************************
         * \brief write - RegistryObject implementation. Field is write down to provided
         *               dictionary(aimed for automatic writing file->dict->object)
         * \param dict - provided dictionaray for data.
         * *******************************************************************************/
        virtual void write(iomanagment::Dictionary &dict)
        {
            using namespace iomanagment;
            
            dict.add(new DictEntry("fieldType",CmpTraits<T>::fieldTypeName()));
            
            DictEntry* internalFieldEntry = new DictEntry("internalField");
            *internalFieldEntry<<*this;
            dict.add(internalFieldEntry);
            
            for (Patch* patch : m_boundaryFields)
            {
                patch->write(dict);
            }
        }
        
        //----------------------------------------------------------//
        // TimeChangeListner implementation           
        //----------------------------------------------------------//
        /** ****************************************************************************
         * \brief timeChanged - listner for time changes
         *******************************************************************************/
        void timeChanged(Scalar time)
        {
            cacheCurrentField();
            
            for(Patch* p : m_boundaryFields)
                p->timeChanged(time);
        }
        
        // ----------------------------------------------------------------//
        //      EXPRESION TEMPLATES HANDLING
        // ----------------------------------------------------------------//
        /// copy constructor from expresion is not allowed --no mesh
        
        /** *****************************************************************************
         * configure asigment operators from expressions
         *******************************************************************************/
        CONFIGURE_ASSIGMENT_OPERATORS_IN_DERIVED_FROM_NUMARRAYBASE( T, GeometricField<T>, FieldBase )
        
        // ----------------------------------------------------------------//
        //      ELEMENT_FILED_BASE ASSIGMENT HANDLING
        // ----------------------------------------------------------------//
        /** ****************************************************************************
         * \function constructor from some discontinued field. Values on element edge
         * are duplicated(not strictly) in discontious field, so this constructor shall
         * be used carefully--> main perpous is to convert disc. field 
         * into some file-writable object. When values are assined, then result value on
         * element edge nodes is undefind, it's some value from 2 connected elements
         * (or more if it's corner)
         *******************************************************************************/
        template<typename EFDerived>
        GeometricField(const ElementFieldBase<T,EFDerived> & ef)
        : FieldBase(ef.mesh().nodesNumber(),"GeometricField"), 
          m_timeRegister(nullptr), m_mesh(ef.mesh())
        {
            for(size_t e=0; e< m_mesh.size(); ++e)
            {
                this->element(e) = ef.element(e);
            }
        }
        
        /** ****************************************************************************
         * \function convert discontinous field into continous result field. Only 
         * assigment is allowed, due to fact, that any other asigment would cause duplication
         * of operation on element edges. Should be used carefully(explenation in above 
         * constructor description)
         * \return reference to this object
         *******************************************************************************/
        template<typename OtherDerived>                                                     
        GeometricField<T> & operator = (const ElementFieldBase<T,OtherDerived> & other)
        {
            for(size_t e=0; e<m_mesh.size(); ++e)
            {
                this->element(e) = other.element(e);
            }
            return *this;
        }
        
    };
    
    template<typename T>
    int GeometricField<T>::NCached = 3;
    
    
    /// \class GeometricFieldElementProxy
    /// Class to allow usage of GeometircField the same as
    /// any ElementFieldBase -- some kinde of discontinous field.
    /// GeometricField is special case of discontinous field kinde, 
    /// but coudn't be directyl derived from ElementFieldBase,
    /// because of it's nature of data storage and defined for that type
    /// mathematical operators. ElementField implements slightly diffrent 
    /// apporach, because there data can be diffrent on the same edge in 
    /// 2 connected edges, so the mathematical simple operators do diffrent
    /// work.
    template<typename T>
    class GeometricFieldElementProxy;
    
    template<typename T>
    struct ElementFieldTraits<T,GeometricFieldElementProxy<T> >
    {
        typedef numArrayIndexMapped< const numArray<T> > return_type;
    };
    
    template<typename T>
    class GeometricFieldElementProxy : public ElementFieldBase<T, GeometricFieldElementProxy<T> >
    {
        const GeometricField<T> & m_gField;
        typedef ElementFieldBase<T, GeometricFieldElementProxy<T> > baseType;
        typedef typename ElementFieldTraits<T,GeometricFieldElementProxy<T>>::return_type return_type;
    public:
        GeometricFieldElementProxy(const GeometricField<T> & gField)
        : baseType(gField.mesh()),m_gField(gField)
        {
        }
        return_type element(size_t e) const { return m_gField.element(e);}
    };
} //field
} //SEM




#endif // GEOMETRICFIELD_H
