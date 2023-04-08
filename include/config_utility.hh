#ifndef CONFIGUTILITY_HH
#define CONFIGUTILITY_HH 1

#include <string>

class config_utility{
    public:
        static config_utility *GetInstance();
        void Initialize(const std::string &configFileName);

        // actual members. Made public for ease of access
        std::string processName;
        std::string nucleusName;
        std::string wfType;
        std::string potentialType;
        std::string potentialRepository;
        std::string potentialFileName;
        std::string wfFileNameSeed;
        std::string outputDirectory;
        double zParent;
        double aParent;
        double maximumRadius;
        double minimumRadius;
        double minimumEnergy;
        double maximumEnergy;
        int nEnergyPoints;
        int nRadialPoints;
        int kBounds[2];
        int maxPrincipalQN;
        int minPrincipalQN;
        bool applyScreening;
        bool writeWF;

        double nuclearRadius;

    private:
        config_utility();
        virtual ~config_utility();
        static config_utility *fInstance;
};

#endif