from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any


class ProviderError(RuntimeError):
    pass


class StructureProvider(ABC):
    name: str
    provider_type: str

    @abstractmethod
    def search_by_elements(
        self,
        elements: list[str],
        max_elements: int,
        include_subsystems: bool,
        filters: dict[str, Any],
    ) -> list[Any]:
        raise NotImplementedError

    @abstractmethod
    def fetch_structure(self, entry_id: str) -> Any:
        raise NotImplementedError


class StubProvider(StructureProvider):
    def __init__(self, name: str, provider_type: str, message: str) -> None:
        self.name = name
        self.provider_type = provider_type
        self.message = message

    def search_by_elements(
        self,
        elements: list[str],
        max_elements: int,
        include_subsystems: bool,
        filters: dict[str, Any],
    ) -> list[Any]:
        raise ProviderError(self.message)

    def fetch_structure(self, entry_id: str) -> Any:
        raise ProviderError(self.message)


class ProviderRegistry:
    def __init__(self) -> None:
        self._providers: dict[str, StructureProvider] = {}

    def register(self, provider: StructureProvider) -> None:
        self._providers[provider.name] = provider

    def get(self, name: str) -> StructureProvider:
        if name not in self._providers:
            raise KeyError(f"provider not registered: {name}")
        return self._providers[name]

    def names(self) -> list[str]:
        return sorted(self._providers)


def default_provider_registry() -> ProviderRegistry:
    registry = ProviderRegistry()
    registry.register(StubProvider("COD", "public_api", "COD provider is scaffolded but not implemented in this build."))
    registry.register(StubProvider("AFLOW", "public_api", "AFLOW provider is scaffolded but not implemented in this build."))
    registry.register(StubProvider("OQMD", "public_api", "OQMD provider is scaffolded but not implemented in this build."))
    registry.register(StubProvider("NOMAD_OPTIMADE", "public_api", "NOMAD/OPTIMADE provider is scaffolded but not implemented in this build."))
    registry.register(StubProvider("ICSD_manual", "restricted_manual", "ICSD is manual-import only. Ophi will not scrape restricted databases."))
    registry.register(StubProvider("ICSD_api", "restricted_api", "ICSD API mode requires explicit legal API credentials."))
    return registry
